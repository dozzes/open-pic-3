#include "core/field_solver.h"
#include "core/cell_stencil.h"
#include "io/io_utilities.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>

#include "pocketfft_hdronly.h"

namespace PIC {

// Faraday's Law predictor step: B(n) -> B(n+1/2)
//   dB/dt = -c * curl(E),  discretised as  B -= c*tau/2/h * curl(E)
//
// B components live on grid edges (Yee stencil); E components on faces.
// ctau_2h = c * tau/2 / h absorbs all scalar prefactors.
//
// With isotropic_curl_enabled each first difference is smoothed across the
// two transverse axes with weights (1-4g, g, g, g, g), g = 1/24 -- the value
// that cancels the leading anisotropic O(h^2*k^4) term of the discrete
// curl-curl operator when the same smoothing is applied to BOTH the Faraday
// curl(E) and the Ampere curl(B). div(curl)=0 is preserved exactly: the
// smoothing factors are scalar operators that commute with the differences.
void calc_magnetic_field_half_time(Grid& grid, const Config::Parameters& params)
{
    print_tm("calc_magnetic_field_half_time");

    const double  ctau_2h = params.ctau_2h();
    const index_t kd      = 1;
    const index_t to_kx   = grid.size_x() - kd;
    const index_t to_ky   = grid.size_y() - kd;
    const index_t to_kz   = grid.size_z() - kd;

    if (params.isotropic_curl_enabled) {
        const double g  = 1.0 / 24.0;
        const double w0 = 1.0 - 4.0 * g;

        // Backward first differences of E components; calling them at shifted
        // indices moves the whole difference to a neighbouring line for the
        // transverse smoothing.
        const auto dx_Ey = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.y - grid(i - 1, j, k).E.y; };
        const auto dx_Ez = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.z - grid(i - 1, j, k).E.z; };
        const auto dy_Ex = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.x - grid(i, j - 1, k).E.x; };
        const auto dy_Ez = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.z - grid(i, j - 1, k).E.z; };
        const auto dz_Ex = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.x - grid(i, j, k - 1).E.x; };
        const auto dz_Ey = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).E.y - grid(i, j, k - 1).E.y; };

#pragma omp parallel for collapse(3)
        for (index_t kx = kd; kx < to_kx; ++kx)
            for (index_t ky = kd; ky < to_ky; ++ky)
                for (index_t kz = kd; kz < to_kz; ++kz) {
                    Cell& cell = grid(kx, ky, kz);

                    // d/dy smoothed over x,z; d/dz smoothed over x,y; etc.
                    const double DyEz = w0 * dy_Ez(kx, ky, kz)
                                        + g * (dy_Ez(kx - 1, ky, kz) + dy_Ez(kx + 1, ky, kz)
                                               + dy_Ez(kx, ky, kz - 1) + dy_Ez(kx, ky, kz + 1));
                    const double DzEy = w0 * dz_Ey(kx, ky, kz)
                                        + g * (dz_Ey(kx - 1, ky, kz) + dz_Ey(kx + 1, ky, kz)
                                               + dz_Ey(kx, ky - 1, kz) + dz_Ey(kx, ky + 1, kz));
                    const double DzEx = w0 * dz_Ex(kx, ky, kz)
                                        + g * (dz_Ex(kx - 1, ky, kz) + dz_Ex(kx + 1, ky, kz)
                                               + dz_Ex(kx, ky - 1, kz) + dz_Ex(kx, ky + 1, kz));
                    const double DxEz = w0 * dx_Ez(kx, ky, kz)
                                        + g * (dx_Ez(kx, ky - 1, kz) + dx_Ez(kx, ky + 1, kz)
                                               + dx_Ez(kx, ky, kz - 1) + dx_Ez(kx, ky, kz + 1));
                    const double DxEy = w0 * dx_Ey(kx, ky, kz)
                                        + g * (dx_Ey(kx, ky - 1, kz) + dx_Ey(kx, ky + 1, kz)
                                               + dx_Ey(kx, ky, kz - 1) + dx_Ey(kx, ky, kz + 1));
                    const double DyEx = w0 * dy_Ex(kx, ky, kz)
                                        + g * (dy_Ex(kx - 1, ky, kz) + dy_Ex(kx + 1, ky, kz)
                                               + dy_Ex(kx, ky, kz - 1) + dy_Ex(kx, ky, kz + 1));

                    cell.B.x -= ctau_2h * (DyEz - DzEy); // Bx at Edge X (i+1/2, j, k)
                    cell.B.y -= ctau_2h * (DzEx - DxEz); // By at Edge Y (i, j+1/2, k)
                    cell.B.z -= ctau_2h * (DxEy - DyEx); // Bz at Edge Z (i, j, k+1/2)
                }
        return;
    }

#pragma omp parallel for collapse(3)
    for (index_t kx = kd; kx < to_kx; ++kx)
        for (index_t ky = kd; ky < to_ky; ++ky)
            for (index_t kz = kd; kz < to_kz; ++kz) {
                const CellStencil s {grid, kx, ky, kz};
                Cell&             cell = s.c();

                // Bx at Edge X (i+1/2, j, k)
                cell.B.x -= ctau_2h * ((cell.E.z - s.ym().E.z) - (cell.E.y - s.zm().E.y));

                // By at Edge Y (i, j+1/2, k)
                cell.B.y -= ctau_2h * ((cell.E.x - s.zm().E.x) - (cell.E.z - s.xm().E.z));

                // Bz at Edge Z (i, j, k+1/2)
                cell.B.z -= ctau_2h * ((cell.E.y - s.xm().E.y) - (cell.E.x - s.ym().E.x));
            }
}

// PSTD (Pseudospectral Time-Domain) Faraday step: B(n) -> B(n+1/2)
//
// Replaces the FDTD finite-difference curl with a spectral curl that uses
// the exact wave-number k instead of the FDTD approximation 2·sin(k·h/2)/h.
// This eliminates the sinc² spatial dispersion error visible with cold electrons.
// Time integration is still leapfrog; the dispersion relation is
//   sin(ω·dt/2) = c·|k|·dt/2  (not ω = c·|k| exactly).
//
// Staggered-grid phase factors:
//   On the Yee C-grid, E lives on faces (j+1/2)*h and B on edges j*h.
//   The spectral derivative of f (sampled at (j+1/2)*h) evaluated at j*h is:
//     DFT[j] -> i*ky * exp(-i*ky*h/2) * F(m)
//   The exp(-i*ky*h/2) shifts the output backward by h/2 to the edge position.
//   Equivalently: ik_stag = i*k*exp(-i*k*h/2) = cd(k*sin(k*h/2), k*cos(k*h/2)).
//
//   Verification: FDTD backward diff (Ez[j]-Ez[j-1])/h in Fourier space gives
//     exp(-i*ky*h/2) * 2i*sin(ky*h/2)/h. PSTD replaces 2*sin(k*h/2)/h -> k
//     (exact), keeping the same exp(-i*k*h/2) phase. At kh->0: ik_stag -> ik.
//     At kh=pi (Nyquist): ik_stag = cd(pi/h, 0), same sign as FDTD.
//
// CFL: determined by the fastest physical wave (Alfvén/magnetosonic), not c.
//   τ·V_fast/h ≤ 2/(√3·π) ≈ 0.367 — enforced in the Lua script via PSTD_tau.
//
// Periodic-boundary limitation:
//   FFT assumes periodic fields.  With absorbing/buffer boundaries the field
//   wraps around spectrally.  set_boundary_conditions() after the step reduces
//   but does not eliminate wrap-around artefacts at open boundaries.
void calc_magnetic_field_pstd_half_time(Grid& grid, const Config::Parameters& params, SolverWorkspace& ws)
{
    print_tm("calc_magnetic_field_pstd_half_time");

    using cd = std::complex<double>;
    using namespace pocketfft;

    const double  ctau_2 = params.ctau_2(); // c * tau/2
    const index_t Nx     = grid.size_x();
    const index_t Ny     = grid.size_y();
    const index_t Nz     = grid.size_z();
    const double  h      = grid.step();
    const size_t  N      = Nx * Ny * Nz;

    // Flat strides for complex [Nx][Ny][Nz] arrays.
    const stride_t sc    = {(ptrdiff_t)(Ny * Nz * sizeof(cd)), (ptrdiff_t)(Nz * sizeof(cd)), (ptrdiff_t)(sizeof(cd))};
    const shape_t  shape = {Nx, Ny, Nz};
    const shape_t  axes  = {0, 1, 2};

    // Complex-to-complex FFT keeps the spectral layout explicit for all axes.
    // PSTD is an experimental solver, so clarity and verification are preferred
    // over the smaller half-complex r2c/c2r storage path.
    // Buffers live in the caller-owned workspace (capacity persists between
    // steps, same allocation behavior as the former function-local statics).
    std::vector<cd>& Ex_r  = ws.pstd_Ex_r;
    std::vector<cd>& Ey_r  = ws.pstd_Ey_r;
    std::vector<cd>& Ez_r  = ws.pstd_Ez_r;
    std::vector<cd>& Ex_c  = ws.pstd_Ex_c;
    std::vector<cd>& Ey_c  = ws.pstd_Ey_c;
    std::vector<cd>& Ez_c  = ws.pstd_Ez_c;
    std::vector<cd>& dBx_c = ws.pstd_dBx_c;
    std::vector<cd>& dBy_c = ws.pstd_dBy_c;
    std::vector<cd>& dBz_c = ws.pstd_dBz_c;
    std::vector<cd>& dBx_r = ws.pstd_dBx_r;
    std::vector<cd>& dBy_r = ws.pstd_dBy_r;
    std::vector<cd>& dBz_r = ws.pstd_dBz_r;
    Ex_r.resize(N);
    Ey_r.resize(N);
    Ez_r.resize(N);
    Ex_c.resize(N);
    Ey_c.resize(N);
    Ez_c.resize(N);
    dBx_c.resize(N);
    dBy_c.resize(N);
    dBz_c.resize(N);
    dBx_r.resize(N);
    dBy_r.resize(N);
    dBz_r.resize(N);

    // Each (i,j,k) writes only its own flat idx across all three loops below:
    // race-free the same way normalize_NP/set_grid_UP are.
#pragma omp parallel for collapse(3)
    for (index_t i = 0; i < Nx; ++i)
        for (index_t j = 0; j < Ny; ++j)
            for (index_t k = 0; k < Nz; ++k) {
                const size_t idx = i * Ny * Nz + j * Nz + k;
                Ex_r[idx]        = cd(grid(i, j, k).E.x, 0.0);
                Ey_r[idx]        = cd(grid(i, j, k).E.y, 0.0);
                Ez_r[idx]        = cd(grid(i, j, k).E.z, 0.0);
            }

    c2c(shape, sc, sc, axes, true, Ex_r.data(), Ex_c.data(), 1.0);
    c2c(shape, sc, sc, axes, true, Ey_r.data(), Ey_c.data(), 1.0);
    c2c(shape, sc, sc, axes, true, Ez_r.data(), Ez_c.data(), 1.0);

    // Spectral curl with staggered-grid phase factors. This is the FDTD
    // backward-difference phase with the finite-difference wavenumber replaced
    // by the exact spectral wavenumber.
    const double dkx = 2.0 * M_PI / (Nx * h);
    const double dky = 2.0 * M_PI / (Ny * h);
    const double dkz = 2.0 * M_PI / (Nz * h);
    const double h2  = 0.5 * h;

    // Not collapse(3): kx/ikx_s (and ky/iky_s below) are hoisted per outer-loop
    // index, so the loop nest isn't perfectly nested as OpenMP collapse requires.
    // Each mi is independent (writes only its own idx range), so parallelizing
    // just the outer loop is race-free and cheap enough per-thread.
#pragma omp parallel for
    for (index_t mi = 0; mi < Nx; ++mi) {
        const double kx = (mi <= Nx / 2) ? dkx * mi : dkx * (static_cast<long>(mi) - static_cast<long>(Nx));
        const cd     ikx_s(kx * std::sin(kx * h2), kx * std::cos(kx * h2));
        for (index_t mj = 0; mj < Ny; ++mj) {
            const double ky = (mj <= Ny / 2) ? dky * mj : dky * (static_cast<long>(mj) - static_cast<long>(Ny));
            const cd     iky_s(ky * std::sin(ky * h2), ky * std::cos(ky * h2));
            for (index_t mk = 0; mk < Nz; ++mk) {
                const double kz = (mk <= Nz / 2) ? dkz * mk : dkz * (static_cast<long>(mk) - static_cast<long>(Nz));
                const cd     ikz_s(kz * std::sin(kz * h2), kz * std::cos(kz * h2));
                const size_t idx = mi * Ny * Nz + mj * Nz + mk;
                const cd     ex = Ex_c[idx], ey = Ey_c[idx], ez = Ez_c[idx];

                // Bx at [HFF]: ∂Ez/∂y|_y=int − ∂Ey/∂z|_z=int
                dBx_c[idx] = -ctau_2 * (iky_s * ez - ikz_s * ey);
                // By at [FHF]: ∂Ex/∂z|_z=int − ∂Ez/∂x|_x=int
                dBy_c[idx] = -ctau_2 * (ikz_s * ex - ikx_s * ez);
                // Bz at [FFH]: ∂Ey/∂x|_x=int − ∂Ex/∂y|_y=int
                dBz_c[idx] = -ctau_2 * (ikx_s * ey - iky_s * ex);
            }
        }
    }

    c2c(shape, sc, sc, axes, false, dBx_c.data(), dBx_r.data(), 1.0 / N);
    c2c(shape, sc, sc, axes, false, dBy_c.data(), dBy_r.data(), 1.0 / N);
    c2c(shape, sc, sc, axes, false, dBz_c.data(), dBz_r.data(), 1.0 / N);

    // Apply ΔB to interior cells (skip boundary layer, same as FDTD).
    // Each (ix,iy,iz) writes only its own cell.B: race-free.
    const index_t kd = 1;
#pragma omp parallel for collapse(3)
    for (index_t ix = kd; ix < Nx - kd; ++ix)
        for (index_t iy = kd; iy < Ny - kd; ++iy)
            for (index_t iz = kd; iz < Nz - kd; ++iz) {
                const size_t idx  = ix * Ny * Nz + iy * Nz + iz;
                Cell&        cell = grid(ix, iy, iz);
                cell.B.x += dBx_r[idx].real();
                cell.B.y += dBy_r[idx].real();
                cell.B.z += dBz_r[idx].real();
            }
}

// Electron velocity from Ampere's law: Ue = Ui - curl(B) / (4π·e·n)
//
// Each component lives on the corresponding face of the Yee (C-grid) stencil.
// The face-averaged density ne_f matches that staggered position.
// Cells filled artificially by set_threshold carry no real plasma;
// the diamagnetic drift is zero there (hybrid approximation breaks down).
//
// With isotropic_curl_enabled this curl(B) uses the same transversely-smoothed
// differences (g = 1/24) as the Faraday curl(E) -- both curls must carry the
// smoothing for the anisotropy cancellation to hold (see
// calc_magnetic_field_half_time).
void calc_electrons_velocity(Grid& grid, const Config::Parameters& params)
{
    print_tm("calc_electrons_velocity");

    const double c_4pi_e_h   = params.c_4pi_e_h();
    const double dens_cutoff = params.dens_cutoff;

    const index_t kd    = 1;
    const index_t to_kx = (grid.size_x() - kd);
    const index_t to_ky = (grid.size_y() - kd);
    const index_t to_kz = (grid.size_z() - kd);

    if (params.isotropic_curl_enabled) {
        const double g  = 1.0 / 24.0;
        const double w0 = 1.0 - 4.0 * g;

        // Forward first differences of B components (edge -> face staggering).
        const auto dx_By = [&grid](index_t i, index_t j, index_t k) { return grid(i + 1, j, k).B.y - grid(i, j, k).B.y; };
        const auto dx_Bz = [&grid](index_t i, index_t j, index_t k) { return grid(i + 1, j, k).B.z - grid(i, j, k).B.z; };
        const auto dy_Bx = [&grid](index_t i, index_t j, index_t k) { return grid(i, j + 1, k).B.x - grid(i, j, k).B.x; };
        const auto dy_Bz = [&grid](index_t i, index_t j, index_t k) { return grid(i, j + 1, k).B.z - grid(i, j, k).B.z; };
        const auto dz_Bx = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k + 1).B.x - grid(i, j, k).B.x; };
        const auto dz_By = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k + 1).B.y - grid(i, j, k).B.y; };

#pragma omp parallel for collapse(3)
        for (index_t kx = kd; kx < to_kx; ++kx)
            for (index_t ky = kd; ky < to_ky; ++ky)
                for (index_t kz = kd; kz < to_kz; ++kz) {
                    const CellStencil s {grid, kx, ky, kz};
                    Cell&             cell = s.c();

                    // --- FACE X: UEx at (i, j+1/2, k+1/2) ---
                    {
                        const double ne_f = 0.5 * (cell.NP + s.xm().NP);
                        const double DyBz = w0 * dy_Bz(kx, ky, kz)
                                            + g * (dy_Bz(kx - 1, ky, kz) + dy_Bz(kx + 1, ky, kz)
                                                   + dy_Bz(kx, ky, kz - 1) + dy_Bz(kx, ky, kz + 1));
                        const double DzBy = w0 * dz_By(kx, ky, kz)
                                            + g * (dz_By(kx - 1, ky, kz) + dz_By(kx + 1, ky, kz)
                                                   + dz_By(kx, ky - 1, kz) + dz_By(kx, ky + 1, kz));
                        const double curl_Bx = DyBz - DzBy;
                        cell.UE.x = (ne_f > dens_cutoff) ? cell.UP.x - (c_4pi_e_h / ne_f) * curl_Bx : 0.0;
                    }

                    // --- FACE Y: UEy at (i+1/2, j, k+1/2) ---
                    {
                        const double ne_f = 0.5 * (cell.NP + s.ym().NP);
                        const double DzBx = w0 * dz_Bx(kx, ky, kz)
                                            + g * (dz_Bx(kx - 1, ky, kz) + dz_Bx(kx + 1, ky, kz)
                                                   + dz_Bx(kx, ky - 1, kz) + dz_Bx(kx, ky + 1, kz));
                        const double DxBz = w0 * dx_Bz(kx, ky, kz)
                                            + g * (dx_Bz(kx, ky - 1, kz) + dx_Bz(kx, ky + 1, kz)
                                                   + dx_Bz(kx, ky, kz - 1) + dx_Bz(kx, ky, kz + 1));
                        const double curl_By = DzBx - DxBz;
                        cell.UE.y = (ne_f > dens_cutoff) ? cell.UP.y - (c_4pi_e_h / ne_f) * curl_By : 0.0;
                    }

                    // --- FACE Z: UEz at (i+1/2, j+1/2, k) ---
                    {
                        const double ne_f = 0.5 * (cell.NP + s.zm().NP);
                        const double DxBy = w0 * dx_By(kx, ky, kz)
                                            + g * (dx_By(kx, ky - 1, kz) + dx_By(kx, ky + 1, kz)
                                                   + dx_By(kx, ky, kz - 1) + dx_By(kx, ky, kz + 1));
                        const double DyBx = w0 * dy_Bx(kx, ky, kz)
                                            + g * (dy_Bx(kx - 1, ky, kz) + dy_Bx(kx + 1, ky, kz)
                                                   + dy_Bx(kx, ky, kz - 1) + dy_Bx(kx, ky, kz + 1));
                        const double curl_Bz = DxBy - DyBx;
                        cell.UE.z = (ne_f > dens_cutoff) ? cell.UP.z - (c_4pi_e_h / ne_f) * curl_Bz : 0.0;
                    }
                }
        return;
    }

#pragma omp parallel for collapse(3)
    for (index_t kx = kd; kx < to_kx; ++kx)
        for (index_t ky = kd; ky < to_ky; ++ky)
            for (index_t kz = kd; kz < to_kz; ++kz) {
                const CellStencil s {grid, kx, ky, kz};
                Cell&             cell = s.c();

                // --- FACE X: UEx at (i, j+1/2, k+1/2) ---
                {
                    const double ne_f = 0.5 * (cell.NP + s.xm().NP);
                    const double curl_Bx = (s.yp().B.z - cell.B.z) - (s.zp().B.y - cell.B.y);
                    cell.UE.x = (ne_f > dens_cutoff) ? cell.UP.x - (c_4pi_e_h / ne_f) * curl_Bx : 0.0;
                }

                // --- FACE Y: UEy at (i+1/2, j, k+1/2) ---
                {
                    const double ne_f = 0.5 * (cell.NP + s.ym().NP);
                    const double curl_By = (s.zp().B.x - cell.B.x) - (s.xp().B.z - cell.B.z);
                    cell.UE.y = (ne_f > dens_cutoff) ? cell.UP.y - (c_4pi_e_h / ne_f) * curl_By : 0.0;
                }

                // --- FACE Z: UEz at (i+1/2, j+1/2, k) ---
                {
                    const double ne_f = 0.5 * (cell.NP + s.zm().NP);
                    const double curl_Bz = (s.xp().B.y - cell.B.y) - (s.yp().B.x - cell.B.x);
                    cell.UE.z = (ne_f > dens_cutoff) ? cell.UP.z - (c_4pi_e_h / ne_f) * curl_Bz : 0.0;
                }
            }
}

// Generalised Ohm's law on the Yee (C-grid) stencil:
//   E = (1/c)*(UE × B) - ∇Pe/(e·ne) + η*(c/4π)*curl(B)
// Cells below dens_cutoff have been filled artificially by set_threshold;
// the hybrid approximation breaks down there, so E is set to zero.
void calc_electric_field(Grid& grid, const Config::Parameters& params)
{
    print_tm("calc_electric_field");

    const double c_4pi       = Constants::c_4pi();
    const double e_charge    = Constants::e();
    const double k_B         = Constants::k_B();
    const double inv_c       = 1.0 / Constants::c();
    const double inv_h       = 1.0 / grid.step();
    const double dens_cutoff = params.dens_cutoff;

    const index_t kd    = 1;
    const index_t to_kx = (grid.size_x() - kd);
    const index_t to_ky = (grid.size_y() - kd);
    const index_t to_kz = (grid.size_z() - kd);

    // Same smoothed curl(B) as calc_electrons_velocity when the isotropic
    // curl is on: the resistive term must use the identical operator, or the
    // eta damping would act on a different discrete mode than the one the
    // solver propagates. g = 0 reduces every D to the plain 2-point diff.
    const bool   iso = params.isotropic_curl_enabled;
    const double g   = iso ? 1.0 / 24.0 : 0.0;
    const double w0  = 1.0 - 4.0 * g;

    const auto dx_By = [&grid](index_t i, index_t j, index_t k) { return grid(i + 1, j, k).B.y - grid(i, j, k).B.y; };
    const auto dx_Bz = [&grid](index_t i, index_t j, index_t k) { return grid(i + 1, j, k).B.z - grid(i, j, k).B.z; };
    const auto dy_Bx = [&grid](index_t i, index_t j, index_t k) { return grid(i, j + 1, k).B.x - grid(i, j, k).B.x; };
    const auto dy_Bz = [&grid](index_t i, index_t j, index_t k) { return grid(i, j + 1, k).B.z - grid(i, j, k).B.z; };
    const auto dz_Bx = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k + 1).B.x - grid(i, j, k).B.x; };
    const auto dz_By = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k + 1).B.y - grid(i, j, k).B.y; };

#pragma omp parallel for collapse(3)
    for (index_t kx = kd; kx < to_kx; ++kx)
        for (index_t ky = kd; ky < to_ky; ++ky)
            for (index_t kz = kd; kz < to_kz; ++kz) {
                const CellStencil s {grid, kx, ky, kz};
                Cell&             cell = s.c();

                // --- FACE X: Ex at (i, j+1/2, k+1/2) ---
                {
                    const double ne_x = 0.5 * (cell.NP + s.xm().NP);
                    const double grad_Pe_x = k_B * (cell.NP * cell.Te - s.xm().NP * s.xm().Te) * inv_h;

                    const double By_f = 0.5 * (cell.B.y + s.zp().B.y);
                    const double Bz_f = 0.5 * (cell.B.z + s.yp().B.z);

                    const double UEy_f = 0.25 * (cell.UE.y + s.yp().UE.y + s.xm().UE.y + s.xm_yp().UE.y);
                    const double UEz_f = 0.25 * (cell.UE.z + s.zp().UE.z + s.xm().UE.z + s.xm_zp().UE.z);

                    const double curl_Bx =
                        !iso ? ((s.yp().B.z - cell.B.z) - (s.zp().B.y - cell.B.y)) * inv_h
                             : (w0 * dy_Bz(kx, ky, kz)
                                + g * (dy_Bz(kx - 1, ky, kz) + dy_Bz(kx + 1, ky, kz)
                                       + dy_Bz(kx, ky, kz - 1) + dy_Bz(kx, ky, kz + 1))
                                - w0 * dz_By(kx, ky, kz)
                                - g * (dz_By(kx - 1, ky, kz) + dz_By(kx + 1, ky, kz)
                                       + dz_By(kx, ky - 1, kz) + dz_By(kx, ky + 1, kz)))
                                   * inv_h;
                    const double eta_f   = 0.5 * (cell.eta + s.xm().eta);

                    cell.E.x = (ne_x > dens_cutoff)
                                   ? (UEz_f * By_f - UEy_f * Bz_f) * inv_c - grad_Pe_x / (e_charge * ne_x)
                                         + eta_f * c_4pi * curl_Bx
                                   : 0.0;
                }

                // --- FACE Y: Ey at (i+1/2, j, k+1/2) ---
                {
                    const double ne_y = 0.5 * (cell.NP + s.ym().NP);
                    const double grad_Pe_y = k_B * (cell.NP * cell.Te - s.ym().NP * s.ym().Te) * inv_h;

                    const double Bz_f = 0.5 * (cell.B.z + s.xp().B.z);
                    const double Bx_f = 0.5 * (cell.B.x + s.zp().B.x);

                    const double UEx_f = 0.25 * (cell.UE.x + s.xp().UE.x + s.ym().UE.x + s.xp_ym().UE.x);
                    const double UEz_f = 0.25 * (cell.UE.z + s.zp().UE.z + s.ym().UE.z + s.ym_zp().UE.z);

                    const double curl_By =
                        !iso ? ((s.zp().B.x - cell.B.x) - (s.xp().B.z - cell.B.z)) * inv_h
                             : (w0 * dz_Bx(kx, ky, kz)
                                + g * (dz_Bx(kx - 1, ky, kz) + dz_Bx(kx + 1, ky, kz)
                                       + dz_Bx(kx, ky - 1, kz) + dz_Bx(kx, ky + 1, kz))
                                - w0 * dx_Bz(kx, ky, kz)
                                - g * (dx_Bz(kx, ky - 1, kz) + dx_Bz(kx, ky + 1, kz)
                                       + dx_Bz(kx, ky, kz - 1) + dx_Bz(kx, ky, kz + 1)))
                                   * inv_h;
                    const double eta_f   = 0.5 * (cell.eta + s.ym().eta);

                    cell.E.y = (ne_y > dens_cutoff)
                                   ? (UEx_f * Bz_f - UEz_f * Bx_f) * inv_c - grad_Pe_y / (e_charge * ne_y)
                                         + eta_f * c_4pi * curl_By
                                   : 0.0;
                }

                // --- FACE Z: Ez at (i+1/2, j+1/2, k) ---
                {
                    const double ne_z = 0.5 * (cell.NP + s.zm().NP);
                    const double grad_Pe_z = k_B * (cell.NP * cell.Te - s.zm().NP * s.zm().Te) * inv_h;

                    const double Bx_f = 0.5 * (cell.B.x + s.yp().B.x);
                    const double By_f = 0.5 * (cell.B.y + s.xp().B.y);

                    const double UEx_f = 0.25 * (cell.UE.x + s.xp().UE.x + s.zm().UE.x + s.xp_zm().UE.x);
                    const double UEy_f = 0.25 * (cell.UE.y + s.yp().UE.y + s.zm().UE.y + s.yp_zm().UE.y);

                    const double curl_Bz =
                        !iso ? ((s.xp().B.y - cell.B.y) - (s.yp().B.x - cell.B.x)) * inv_h
                             : (w0 * dx_By(kx, ky, kz)
                                + g * (dx_By(kx, ky - 1, kz) + dx_By(kx, ky + 1, kz)
                                       + dx_By(kx, ky, kz - 1) + dx_By(kx, ky, kz + 1))
                                - w0 * dy_Bx(kx, ky, kz)
                                - g * (dy_Bx(kx - 1, ky, kz) + dy_Bx(kx + 1, ky, kz)
                                       + dy_Bx(kx, ky, kz - 1) + dy_Bx(kx, ky, kz + 1)))
                                   * inv_h;
                    const double eta_f   = 0.5 * (cell.eta + s.zm().eta);

                    cell.E.z = (ne_z > dens_cutoff)
                                   ? (UEy_f * Bx_f - UEx_f * By_f) * inv_c - grad_Pe_z / (e_charge * ne_z)
                                         + eta_f * c_4pi * curl_Bz
                                   : 0.0;
                }
            }

}

} // namespace PIC
