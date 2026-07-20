#pragma once

#include "io/io_utilities.h"
#include "grid/grid.h"
#include "config/config.h"
#include <vector>
#include <algorithm>

namespace PIC {

namespace detail {

// Implicit tridiagonal sweep (Thomas algorithm) for one grid line of the heat
// equation: (I - coeff*D[chi_face]D) T_new = rhs, chi at faces as harmonic
// mean. Endpoints [first-1] and [last+1] enter the RHS as Dirichlet values
// taken from T_old. Free function (not a HeatSolver member) so unit tests can
// verify it against a dense solve.
inline void run_Thomas(const std::vector<double>& chi,
                       const std::vector<double>& T_old,
                       std::vector<double>&       T_new,
                       int                        first,
                       int                        last,
                       double                     coeff,
                       std::vector<double>&       p,
                       std::vector<double>&       q)
{
    if (last < first)
        return;

    // p/q are caller-owned scratch: run_Thomas fires once per grid line
    // (nx*ny lines per sweep), so allocating here would dominate the sweep.
    const int K = last - first + 1;
    p.resize(K);
    q.resize(K);

    // HARMONIC MEAN for chi at cell faces (ensures flux symmetry).
    auto get_chi_face = [&](int k1, int k2) {
        return 2.0 * (chi[k1] * chi[k2]) / (chi[k1] + chi[k2] + 1e-20);
    };

    {
        const int    k   = first;
        const double ak  = -coeff * get_chi_face(k - 1, k);
        const double ck  = -coeff * get_chi_face(k, k + 1);
        const double bk  = 1.0 - ak - ck;
        const double rhs = T_old[k] - ak * T_old[k - 1];
        p[0]             = -ck / bk;
        q[0]             = rhs / bk;
    }

    for (int k = first + 1; k < last; ++k) {
        const int n  = k - first;
        double    ak = -coeff * get_chi_face(k - 1, k);
        double    ck = -coeff * get_chi_face(k, k + 1);
        double    bk = 1.0 - ak - ck;
        double    m  = bk + ak * p[n - 1];
        p[n]         = -ck / m;
        q[n]         = (T_old[k] - ak * q[n - 1]) / m;
    }

    if (K == 1) {
        const int    k  = first;
        const double ak = -coeff * get_chi_face(k - 1, k);
        const double ck = -coeff * get_chi_face(k, k + 1);
        const double bk = 1.0 - ak - ck;
        T_new[k]        = (T_old[k] - ak * T_old[k - 1] - ck * T_old[k + 1]) / bk;
        return;
    }

    {
        const int k   = last;
        const int n   = K - 1;
        double    ak  = -coeff * get_chi_face(k - 1, k);
        double    ck  = -coeff * get_chi_face(k, k + 1);
        double    bk  = 1.0 - ak - ck;
        double    rhs = T_old[k] - ck * T_old[k + 1];
        T_new[k]      = (rhs - ak * q[n - 1]) / (bk + ak * p[n - 1]);
    }

    for (int k = last - 1; k >= first; --k) {
        const int n = k - first;
        T_new[k]    = p[n] * T_new[k + 1] + q[n];
    }
}

} // namespace detail

class HeatSolver
{
  public:
    // `step` is the global time-step number: solve() runs exactly once per
    // step (starting at 1), so using it for the sweep alternation reproduces
    // the former hidden static call counter -- without state that silently
    // carried over between runs in the same process.
    static void solve(Grid& grid, index_t step)
    {
        print_tm("HeatSolver::solve");

        double chi_base = Config::electron_thermal_conductivity();
        if (chi_base <= 0.0)
            return;

        const double dt    = Config::tau() / 3.0;
        const double h     = grid.step();
        const double coeff = (chi_base * dt) / (h * h);

        const int nx = (int)grid.size_x();
        const int ny = (int)grid.size_y();
        const int nz = (int)grid.size_z();

        // 1. Snapshot of thermal conductivity with vacuum-cell protection.
        std::vector<double> chi_snap((size_t)nx * ny * nz);
        const double        dens_cutoff = Config::dens_cutoff();
#pragma omp parallel for collapse(3) schedule(static)
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    const Cell& c = grid(ix, iy, iz);
                    // No plasma → no thermal conductivity.
                    // Using dens_cutoff (same floor as set_threshold) prevents
                    // heat diffusion into vacuum cells and Te runaway there.
                    if (c.NP < dens_cutoff) {
                        chi_snap[(ix * ny + iy) * nz + iz] = 0.0;
                    } else {
                        double local_chi                   = c.Te / (c.eta + 1e-12);
                        chi_snap[(ix * ny + iy) * nz + iz] = std::clamp(local_chi, 1e-6, 100.0);
                    }
                }
            }
        }

        // Sweep direction alternates each step to suppress directional bias.
        const bool flip = (step % 2 != 0);

        // Axis order also alternates each step.
        if (step % 2 == 0) {
            solve_sweep(grid, chi_snap, coeff, 0, flip);
            solve_sweep(grid, chi_snap, coeff, 1, flip);
            solve_sweep(grid, chi_snap, coeff, 2, flip);
        } else {
            solve_sweep(grid, chi_snap, coeff, 2, flip);
            solve_sweep(grid, chi_snap, coeff, 1, flip);
            solve_sweep(grid, chi_snap, coeff, 0, flip);
        }
    }

  private:
    static void solve_sweep(Grid& grid, const std::vector<double>& chi_snap, double coeff, int dir, bool flip)
    {
        const int nx = (int)grid.size_x();
        const int ny = (int)grid.size_y();
        const int nz = (int)grid.size_z();

        int K_max, I_max, J_max;
        if (dir == 0) {
            K_max = nx;
            I_max = ny;
            J_max = nz;
        } else if (dir == 1) {
            K_max = ny;
            I_max = nx;
            J_max = nz;
        } else {
            K_max = nz;
            I_max = nx;
            J_max = ny;
        }

#pragma omp parallel
        {
            std::vector<double> chi_l(K_max), T_old(K_max), T_new(K_max);
            std::vector<double> p_buf(K_max), q_buf(K_max); // run_Thomas scratch, reused per line
#pragma omp for collapse(2) schedule(static)
            for (int i = 2; i < I_max - 2; ++i) {
                for (int j = 2; j < J_max - 2; ++j) {
                    for (int k = 0; k < K_max; ++k) {
                        int ix, iy, iz;
                        if (dir == 0) {
                            ix = k;
                            iy = i;
                            iz = j;
                        } else if (dir == 1) {
                            ix = i;
                            iy = k;
                            iz = j;
                        } else {
                            ix = i;
                            iy = j;
                            iz = k;
                        }

                        T_old[k] = grid(ix, iy, iz).Te;
                        chi_l[k] = chi_snap[(ix * ny + iy) * nz + iz];
                    }

                    if (flip) {
                        std::reverse(chi_l.begin(), chi_l.end());
                        std::reverse(T_old.begin(), T_old.end());
                    }

                    T_new = T_old;

                    const int first = 2;
                    const int last  = K_max - 3;
                    detail::run_Thomas(chi_l, T_old, T_new, first, last, coeff, p_buf, q_buf);

                    if (flip) {
                        std::reverse(T_new.begin(), T_new.end());
                    }

                    for (int k = 2; k < K_max - 2; ++k) {
                        int ix, iy, iz;
                        if (dir == 0) {
                            ix = k;
                            iy = i;
                            iz = j;
                        } else if (dir == 1) {
                            ix = i;
                            iy = k;
                            iz = j;
                        } else {
                            ix = i;
                            iy = j;
                            iz = k;
                        }
                        grid(ix, iy, iz).Te = std::max(T_new[k], 1e-4);
                    }
                }
            }
        }
    }
};
} // namespace PIC
