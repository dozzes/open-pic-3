#include "config/config.h"
#include "config/constants.h"
#include "config/diagnostics.h"
#include "core/gather_scatter.h"
#include "grid/grid.h"
#include "core/heat_solver.h"
#include "mpi/mpi_runtime.h"
#include "particles/particle_groups.h"
#include "particles/particles.h"
#include "core/simulate.h"
#include "lua/use_lua.h"
#include "util/utility.h"
#include "pocketfft_hdronly.h"

#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void expect_true(const std::string& name, bool condition);

void expect_near(const std::string& name, double actual, double expected, double abs_tol, double rel_tol = 0.0)
{
    const double scale = std::max(std::fabs(expected), 1.0);
    const double tol   = std::max(abs_tol, rel_tol * scale);
    if (std::fabs(actual - expected) > tol) {
        std::cerr << "FAIL " << name << ": actual=" << actual << " expected=" << expected << " tol=" << tol << "\n";
        ++g_failures;
    }
}

void expect_vector_near(const std::string& name,
                        const DblVector&  actual,
                        const DblVector&  expected,
                        double            abs_tol,
                        double            rel_tol = 0.0)
{
    expect_near(name + ".x", actual.x, expected.x, abs_tol, rel_tol);
    expect_near(name + ".y", actual.y, expected.y, abs_tol, rel_tol);
    expect_near(name + ".z", actual.z, expected.z, abs_tol, rel_tol);
}

int reset_single_ion_group()
{
    ParticleGroups groups;
    groups.clear();
    groups.create_group("ion", 1.0, 1.0, PIC::no_diag);
    return ParticleGroups::get_id_by_name("ion");
}

void reset_test_config(double h, double tau)
{
    // Config is a process-wide singleton: rebuild it from defaults so a field
    // forgotten here cannot leak from the previous test.
    PIC::Config::parameters()     = PIC::Config::Parameters();
    PIC::Config::parameters().h   = h;
    PIC::Config::parameters().tau = tau;
}

Grid make_uniform_grid(const DblVector& e, const DblVector& b)
{
    Grid grid(8, 8, 8, PIC::Config::h());
    for (index_t i = 0; i < grid.size_x(); ++i) {
        for (index_t j = 0; j < grid.size_y(); ++j) {
            for (index_t k = 0; k < grid.size_z(); ++k) {
                grid(i, j, k).E = e;
                grid(i, j, k).B = b;
            }
        }
    }
    return grid;
}

void test_Boris_uniform_e()
{
    reset_test_config(1.0, 1.0e-15);
    const int group_id = reset_single_ion_group();

    const DblVector e(0.75, -0.25, 0.5);
    const DblVector b(0.0, 0.0, 0.0);
    const Grid grid = make_uniform_grid(e, b);

    const DblVector v0(1.0e5, -2.0e5, 3.0e5);
    Particle particle(group_id, DblVector(3.25, 3.35, 3.45), v0, 1.0);

    PIC::push_particle_Boris(grid, particle);

    const double qm = PIC::Constants::e_mp();
    const DblVector expected = v0 + PIC::Config::tau() * qm * e;
    expect_vector_near("Boris_uniform_e.velocity", particle.v, expected, 1.0e-9, 1.0e-14);
}

void test_Boris_uniform_b_preserves_speed()
{
    reset_test_config(1.0, 1.0e-12);
    const int group_id = reset_single_ion_group();

    const DblVector e(0.0, 0.0, 0.0);
    const DblVector b(0.0, 0.0, 100.0);
    const Grid grid = make_uniform_grid(e, b);

    Particle particle(group_id, DblVector(3.25, 3.35, 3.45), DblVector(1.0e5, 2.0e5, -0.5e5), 1.0);
    const double speed0 = particle.v.get_len();

    for (int i = 0; i < 1000; ++i) {
        PIC::push_particle_Boris(grid, particle);
    }

    expect_near("Boris_uniform_b.speed", particle.v.get_len(), speed0, 1.0e-7, 1.0e-12);
}

void test_cic_scatter_weight_sum()
{
    reset_test_config(1.0, 1.0e-12);
    const int group_id = reset_single_ion_group();

    const std::vector<DblVector> positions = {
        DblVector(3.5, 3.5, 3.5),
        DblVector(3.0, 3.5, 3.5),
        DblVector(3.23, 3.41, 3.67),
    };

    for (size_t case_idx = 0; case_idx < positions.size(); ++case_idx) {
        PIC::DensityGrid density_grid(8, 8, 8, PIC::Config::h());
        const DblVector  velocity(1.2, -0.4, 0.8);
        const double     ni = 2.5;
        const Particle   particle(group_id, positions[case_idx], velocity, ni);

        PIC::scatter_particle_snapshot(particle, density_grid);

        double    np_sum = 0.0;
        DblVector up_sum(0.0, 0.0, 0.0);
        DblVector up_np_sum(0.0, 0.0, 0.0);
        for (size_t i = 0; i < density_grid.raw_size(); ++i) {
            const Density& cell = density_grid.raw_get(i);
            np_sum += cell.NP;
            up_sum += cell.UP;
            up_np_sum += cell.UP_NP;
        }

        const std::string prefix = "cic_case_" + std::to_string(case_idx);
        expect_near(prefix + ".NP_sum", np_sum, ni, 1.0e-12, 1.0e-13);
        expect_vector_near(prefix + ".UP_sum", up_sum, ni * velocity, 1.0e-12, 1.0e-13);
        expect_vector_near(prefix + ".UP_NP_sum", up_np_sum, DblVector(ni, ni, ni), 1.0e-12, 1.0e-13);
    }
}

void test_ngp_staggered_gather_and_scatter()
{
    reset_test_config(1.0, 1.0e-12);
    PIC::Config::parameters().scatter_alg = PIC::ScatterAlg::NGP;
    const int group_id = reset_single_ion_group();

    Grid grid(8, 8, 8, 1.0);
    for (index_t i = 0; i < grid.size_x(); ++i)
        for (index_t j = 0; j < grid.size_y(); ++j)
            for (index_t k = 0; k < grid.size_z(); ++k) {
                const double code = 100.0 * i + 10.0 * j + k;
                grid(i, j, k).NP = code;
                grid(i, j, k).E  = DblVector(code, code, code);
                grid(i, j, k).B  = DblVector(code, code, code);
            }

    // The deliberately asymmetric z coordinate makes the nearest integer
    // lattice point differ from the nearest half-integer lattice point.
    const DblVector p(3.23, 3.41, 3.67);
    DblVector       gathered;

    PIC::gather_face(grid, p, &Cell::E, gathered);
    expect_vector_near("ngp_gather.face", gathered, DblVector(333.0, 333.0, 334.0), 0.0);

    PIC::gather_edge(grid, p, &Cell::B, gathered);
    expect_vector_near("ngp_gather.edge", gathered, DblVector(334.0, 334.0, 333.0), 0.0);

    expect_near("ngp_gather.center", PIC::gather_center(grid, p, &Cell::NP), 333.0, 0.0);

    const DblVector velocity(1.2, -0.4, 0.8);
    const double    ni = 2.5;
    const Particle  particle(group_id, p, velocity, ni);
    PIC::DensityGrid density_grid(8, 8, 8, 1.0);
    PIC::scatter_particle_snapshot(particle, density_grid);

    expect_near("ngp_scatter.NP", density_grid(3, 3, 3).NP, ni, 0.0);
    expect_near("ngp_scatter.UPx", density_grid(3, 3, 3).UP.x, ni * velocity.x, 0.0);
    expect_near("ngp_scatter.UPy", density_grid(3, 3, 3).UP.y, ni * velocity.y, 0.0);
    expect_near("ngp_scatter.UPz", density_grid(3, 3, 4).UP.z, ni * velocity.z, 0.0);
    expect_near("ngp_scatter.UP_NPx", density_grid(3, 3, 3).UP_NP.x, ni, 0.0);
    expect_near("ngp_scatter.UP_NPy", density_grid(3, 3, 3).UP_NP.y, ni, 0.0);
    expect_near("ngp_scatter.UP_NPz", density_grid(3, 3, 4).UP_NP.z, ni, 0.0);

    double    np_sum = 0.0;
    DblVector up_sum(0.0, 0.0, 0.0);
    DblVector up_np_sum(0.0, 0.0, 0.0);
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        const Density& node = density_grid.raw_get(i);
        np_sum += node.NP;
        up_sum += node.UP;
        up_np_sum += node.UP_NP;
    }
    expect_near("ngp_scatter.NP_sum", np_sum, ni, 0.0);
    expect_vector_near("ngp_scatter.UP_sum", up_sum, ni * velocity, 0.0);
    expect_vector_near("ngp_scatter.UP_NP_sum", up_np_sum, DblVector(ni, ni, ni), 0.0);
}

void test_tsc_staggered_gather_and_scatter()
{
    reset_test_config(1.0, 1.0e-12);
    PIC::Config::parameters().scatter_alg = PIC::ScatterAlg::TSC;
    const int group_id = reset_single_ion_group();

    const DblVector field_value(2.0, -3.0, 4.0);
    Grid            grid(8, 8, 8, 1.0);
    for (index_t i = 0; i < grid.size_x(); ++i)
        for (index_t j = 0; j < grid.size_y(); ++j)
            for (index_t k = 0; k < grid.size_z(); ++k) {
                grid(i, j, k).NP = 7.0;
                grid(i, j, k).E  = field_value;
                grid(i, j, k).B  = field_value;
            }

    const DblVector p(3.23, 3.41, 3.67);
    DblVector       gathered;
    PIC::gather_face(grid, p, &Cell::E, gathered);
    expect_vector_near("tsc_gather.face", gathered, field_value, 1.0e-12, 1.0e-13);
    PIC::gather_edge(grid, p, &Cell::B, gathered);
    expect_vector_near("tsc_gather.edge", gathered, field_value, 1.0e-12, 1.0e-13);
    expect_near("tsc_gather.center", PIC::gather_center(grid, p, &Cell::NP), 7.0, 1.0e-12, 1.0e-13);

    const DblVector velocity(1.2, -0.4, 0.8);
    const double    ni = 2.5;
    const Particle  particle(group_id, p, velocity, ni);
    PIC::DensityGrid density_grid(8, 8, 8, 1.0);
    PIC::scatter_particle_snapshot(particle, density_grid);

    double    np_sum = 0.0;
    DblVector up_sum(0.0, 0.0, 0.0);
    DblVector up_np_sum(0.0, 0.0, 0.0);
    int       np_nodes = 0;
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        const Density& node = density_grid.raw_get(i);
        np_sum += node.NP;
        up_sum += node.UP;
        up_np_sum += node.UP_NP;
        if (node.NP > 0.0)
            ++np_nodes;
    }
    expect_near("tsc_scatter.NP_sum", np_sum, ni, 1.0e-12, 1.0e-13);
    expect_vector_near("tsc_scatter.UP_sum", up_sum, ni * velocity, 1.0e-12, 1.0e-13);
    expect_vector_near("tsc_scatter.UP_NP_sum", up_np_sum, DblVector(ni, ni, ni), 1.0e-12, 1.0e-13);
    expect_true("tsc_scatter.27_density_nodes", np_nodes == 27);

    // Discrete adjoint check: depositing a particle moment and taking its dot
    // product with an arbitrary grid field must equal gathering that field to
    // the particle and multiplying by the same particle moment.
    Grid adjoint_grid(8, 8, 8, 1.0);
    for (index_t i = 0; i < adjoint_grid.size_x(); ++i)
        for (index_t j = 0; j < adjoint_grid.size_y(); ++j)
            for (index_t k = 0; k < adjoint_grid.size_z(); ++k) {
                const double code = 0.7 * i - 0.2 * j + 0.13 * k + 0.01 * i * j;
                adjoint_grid(i, j, k).NP = code;
                adjoint_grid(i, j, k).E = DblVector(code + 1.0, 2.0 * code - 0.5, -0.3 * code + 0.2);
            }
    DblVector gathered_e;
    PIC::gather_face(adjoint_grid, p, &Cell::E, gathered_e);
    const double gathered_np = PIC::gather_center(adjoint_grid, p, &Cell::NP);
    double       np_dot = 0.0;
    DblVector    up_dot(0.0, 0.0, 0.0);
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        const Density& deposited = density_grid.raw_get(i);
        const Cell&    field = adjoint_grid.raw_get(i);
        np_dot += deposited.NP * field.NP;
        up_dot.x += deposited.UP.x * field.E.x;
        up_dot.y += deposited.UP.y * field.E.y;
        up_dot.z += deposited.UP.z * field.E.z;
    }
    expect_near("tsc_adjoint.NP", np_dot, ni * gathered_np, 1.0e-12, 1.0e-13);
    expect_vector_near("tsc_adjoint.UP", up_dot,
                       DblVector(ni * velocity.x * gathered_e.x,
                                 ni * velocity.y * gathered_e.y,
                                 ni * velocity.z * gathered_e.z),
                       1.0e-12, 1.0e-13);
}

void test_pstd_plane_wave_Faraday_derivative()
{
    const double h   = 1.0;
    const double tau = 1.0e-13;
    reset_test_config(h, tau);

    constexpr index_t n = 8;
    Grid              grid(n, n, n, h);

    const int    mode = 1;
    const double ky   = 2.0 * PIC::Constants::pi() * mode / (n * h);
    expect_near("pstd_config.ctau_2", PIC::Config::ctau_2(), PIC::Constants::c() * tau * 0.5, 1.0e-18, 1.0e-14);

    for (index_t i = 0; i < n; ++i) {
        for (index_t j = 0; j < n; ++j) {
            for (index_t k = 0; k < n; ++k) {
                grid(i, j, k).E = DblVector(0.0, 0.0, std::sin(ky * (j + 0.5) * h));
                grid(i, j, k).B = DblVector(0.0, 0.0, 0.0);
            }
        }
    }

    PIC::SolverWorkspace pstd_workspace;
    PIC::calc_magnetic_field_pstd_half_time(grid, PIC::Config::parameters(), pstd_workspace);

    const double ctau_2 = PIC::Constants::c() * tau * 0.5;
    for (index_t i = 1; i < n - 1; ++i) {
        for (index_t j = 1; j < n - 1; ++j) {
            for (index_t k = 1; k < n - 1; ++k) {
                const double expected_bx = -ctau_2 * ky * std::cos(ky * j * h);
                const std::string prefix = "pstd_plane_wave.Bx_" + std::to_string(i) + "_" + std::to_string(j) + "_"
                                           + std::to_string(k);
                expect_near(prefix, grid(i, j, k).B.x, expected_bx, 1.0e-12, 1.0e-10);
                expect_near(prefix + ".By", grid(i, j, k).B.y, 0.0, 1.0e-12);
                expect_near(prefix + ".Bz", grid(i, j, k).B.z, 0.0, 1.0e-12);
            }
        }
    }
}

// Numerical-dispersion symbol of the FDTD Faraday step, measured exactly.
//
// The solver is local, so feeding it a single-component plane wave
// E ~ sin(k·r) (sampled at that component's staggered position) must return
// dB in which every first difference carries the Yee dispersion factor
// 2*sin(k_a*h/2)/h in place of the exact wavenumber k_a, and -- with the
// isotropized curl -- additionally the transverse smoothing factor
// s_a = 1 - 4g*[sin^2(k_b*h/2) + sin^2(k_c*h/2)], g = 1/24. Machine-precision
// agreement pins down the full dispersion symbol of the implementation: a
// wrong sign, index, or smoothing axis fails the test. The three cases below
// exercise all six first differences of the curl.
//
// From the verified symbol the test then derives the whistler-branch
// (omega ~ K^2) anisotropy numbers used in the m=4 front-mode analysis:
// the axis-vs-body-diagonal spread of K^2/k^2 is (kh)^2/18 for the plain
// curl, is cancelled to O((kh)^4) by the isotropized curl, and grows to the
// irreducible face-diagonal/axis ratio 2*sin^2(pi/(2*sqrt(2))) ~ 1.61 at
// the grid scale |k| = pi/h.
void test_fdtd_curl_dispersion_symbol()
{
    const double h = 1.0;
    // tau such that ctau_2h = c*tau/(2h) = 1: dB then equals the raw
    // (smoothed) curl of E and the symbol is read off directly.
    reset_test_config(h, 2.0 * h / PIC::Constants::c());

    constexpr index_t n  = 24;
    const double      k1 = 0.6 / h; // kh = 0.6 and 1.1: the front-relevant
    const double      k2 = 1.1 / h; // scales where the m=4 mode lives

    // One E component per case, oblique wave in its two curl-relevant axes.
    struct WaveCase {
        int    comp;       // 0 = Ex, 1 = Ey, 2 = Ez
        double kx, ky, kz;
    };
    const WaveCase wave_cases[] = {
        {2, k1, k2, 0.0}, // Ez: exercises dx_Ez, dy_Ez
        {0, 0.0, k1, k2}, // Ex: exercises dy_Ex, dz_Ex
        {1, k2, 0.0, k1}, // Ey: exercises dx_Ey, dz_Ey
    };

    for (int iso = 0; iso <= 1; ++iso) {
        PIC::Config::parameters().isotropic_curl_enabled = (iso == 1);
        const double g = (iso == 1) ? 1.0 / 24.0 : 0.0;
        // Transverse smoothing factor of a first difference whose two
        // transverse wavenumbers are kb, kc.
        const auto s_of = [g, h](double kb, double kc) {
            const double sb = std::sin(0.5 * kb * h);
            const double sc = std::sin(0.5 * kc * h);
            return 1.0 - 4.0 * g * (sb * sb + sc * sc);
        };
        const auto yee = [h](double ka) { return 2.0 * std::sin(0.5 * ka * h); };

        for (const WaveCase& wc : wave_cases) {
            Grid grid(n, n, n, h);
            for (index_t i = 0; i < n; ++i)
                for (index_t j = 0; j < n; ++j)
                    for (index_t k = 0; k < n; ++k) {
                        // Staggered positions: Ex at (i, j+1/2, k+1/2),
                        // Ey at (i+1/2, j, k+1/2), Ez at (i+1/2, j+1/2, k).
                        double x = i * h, y = j * h, z = k * h;
                        if (wc.comp == 0) { y += 0.5 * h; z += 0.5 * h; }
                        if (wc.comp == 1) { x += 0.5 * h; z += 0.5 * h; }
                        if (wc.comp == 2) { x += 0.5 * h; y += 0.5 * h; }
                        const double e = std::sin(wc.kx * x + wc.ky * y + wc.kz * z);
                        Cell& cell = grid(i, j, k);
                        cell.E = DblVector(wc.comp == 0 ? e : 0.0,
                                           wc.comp == 1 ? e : 0.0,
                                           wc.comp == 2 ? e : 0.0);
                        cell.B = DblVector(0.0, 0.0, 0.0);
                    }

            PIC::calc_magnetic_field_half_time(grid, PIC::Config::parameters());

            const std::string tag = "fdtd_dispersion.iso" + std::to_string(iso)
                                    + ".comp" + std::to_string(wc.comp) + ".";
            for (index_t i = 1; i + 1 < n; ++i)
                for (index_t j = 1; j + 1 < n; ++j)
                    for (index_t k = 1; k + 1 < n; ++k) {
                        // dB = -curl(E) per component; each difference lands
                        // at the B position (half-cell back along its axis)
                        // and carries yee(k_a) and its smoothing factor.
                        double ebx = 0.0, eby = 0.0, ebz = 0.0;
                        if (wc.comp == 2) {
                            ebx = -s_of(wc.kx, 0.0) * yee(wc.ky)
                                  * std::cos(wc.kx * (i + 0.5) * h + wc.ky * j * h);
                            eby = s_of(wc.ky, 0.0) * yee(wc.kx)
                                  * std::cos(wc.kx * i * h + wc.ky * (j + 0.5) * h);
                        } else if (wc.comp == 0) {
                            eby = -s_of(wc.ky, 0.0) * yee(wc.kz)
                                  * std::cos(wc.ky * (j + 0.5) * h + wc.kz * k * h);
                            ebz = s_of(wc.kz, 0.0) * yee(wc.ky)
                                  * std::cos(wc.ky * j * h + wc.kz * (k + 0.5) * h);
                        } else {
                            ebx = s_of(wc.kx, 0.0) * yee(wc.kz)
                                  * std::cos(wc.kx * (i + 0.5) * h + wc.kz * k * h);
                            ebz = -s_of(wc.kz, 0.0) * yee(wc.kx)
                                  * std::cos(wc.kx * i * h + wc.kz * (k + 0.5) * h);
                        }
                        const Cell& cell = grid(i, j, k);
                        const std::string at = std::to_string(i) + "_" + std::to_string(j)
                                               + "_" + std::to_string(k);
                        expect_near(tag + "Bx_" + at, cell.B.x, ebx, 1.0e-12, 1.0e-12);
                        expect_near(tag + "By_" + at, cell.B.y, eby, 1.0e-12, 1.0e-12);
                        expect_near(tag + "Bz_" + at, cell.B.z, ebz, 1.0e-12, 1.0e-12);
                    }
        }
    }
    PIC::Config::parameters().isotropic_curl_enabled = false;

    // Whistler-branch symbol K^2 = sum_a [s_a * (2/h) sin(k_a h/2)]^2 built
    // from the factors verified above; R = K^2/k^2 is the dispersion error.
    const auto K2_over_k2 = [h](double kx, double ky, double kz, double g) {
        const auto d = [h](double ka) { return 2.0 / h * std::sin(0.5 * ka * h); };
        const auto s = [g, h](double kb, double kc) {
            const double sb = std::sin(0.5 * kb * h);
            const double sc = std::sin(0.5 * kc * h);
            return 1.0 - 4.0 * g * (sb * sb + sc * sc);
        };
        const double Dx = s(ky, kz) * d(kx);
        const double Dy = s(kz, kx) * d(ky);
        const double Dz = s(kx, ky) * d(kz);
        return (Dx * Dx + Dy * Dy + Dz * Dz) / (kx * kx + ky * ky + kz * kz);
    };

    // Long-wave anisotropy: body diagonal vs axis at |k|h = 0.3. The error
    // is LARGEST along a grid axis (each component carries the full k):
    // R_axis ~ 1-(kh)^2/12, R_diag ~ 1-(kh)^2/36, so waves along the body
    // diagonal are the fastest -- the classic FDTD result.
    const double kh = 0.3;
    const double kd = kh / std::sqrt(3.0);
    const double spread_plain = K2_over_k2(kd / h, kd / h, kd / h, 0.0)
                                - K2_over_k2(kh / h, 0.0, 0.0, 0.0);
    // Whistler omega ~ K^2, so this spread IS delta-omega/omega = (kh)^2/18.
    expect_near("fdtd_dispersion.spread_plain_kh0.3", spread_plain, kh * kh / 18.0,
                0.1 * kh * kh / 18.0);
    const double spread_iso = K2_over_k2(kd / h, kd / h, kd / h, 1.0 / 24.0)
                              - K2_over_k2(kh / h, 0.0, 0.0, 1.0 / 24.0);
    expect_true("fdtd_dispersion.iso_cancels_leading_term",
                std::fabs(spread_iso) < 0.05 * std::fabs(spread_plain));

    // Grid-scale anisotropy at |k| = pi/h: face diagonal vs axis. This is the
    // irreducible order-unity ratio no small-kh stencil correction removes.
    const double pi_h  = PIC::Constants::pi() / h;
    const double ratio = K2_over_k2(pi_h / std::sqrt(2.0), pi_h / std::sqrt(2.0), 0.0, 0.0)
                         / K2_over_k2(pi_h, 0.0, 0.0, 0.0);
    const double exact = 2.0 * std::pow(std::sin(PIC::Constants::pi() / (2.0 * std::sqrt(2.0))), 2);
    expect_near("fdtd_dispersion.nyquist_ratio_closed_form", ratio, exact, 1.0e-12, 1.0e-12);
    expect_near("fdtd_dispersion.nyquist_ratio_1.61", ratio, 1.61, 0.01);
}

void test_pocketfft_c2c_roundtrip()
{
    using cd = std::complex<double>;
    using namespace pocketfft;

    const shape_t  shape = {4, 4, 4};
    const stride_t stride = {
        static_cast<ptrdiff_t>(4 * 4 * sizeof(cd)),
        static_cast<ptrdiff_t>(4 * sizeof(cd)),
        static_cast<ptrdiff_t>(sizeof(cd)),
    };
    const shape_t axes = {0, 1, 2};

    std::vector<cd> input(64);
    std::vector<cd> spectrum(64);
    std::vector<cd> output(64);
    for (size_t i = 0; i < input.size(); ++i)
        input[i] = cd(std::sin(0.25 * static_cast<double>(i)), 0.0);

    c2c(shape, stride, stride, axes, true, input.data(), spectrum.data(), 1.0);
    c2c(shape, stride, stride, axes, false, spectrum.data(), output.data(), 1.0 / input.size());

    for (size_t i = 0; i < input.size(); ++i)
        expect_near("pocketfft_c2c_roundtrip_" + std::to_string(i), output[i].real(), input[i].real(), 1.0e-12, 1.0e-12);
}

void expect_true(const std::string& name, bool condition)
{
    if (!condition) {
        std::cerr << "FAIL " << name << "\n";
        ++g_failures;
    }
}

// The unified simulate() driver (serial + MPI merged) calls Mpi::allreduce_*
// and gates I/O on Mpi::is_root() unconditionally. Serial runs stay bit-exact
// only if, without MPI_Init: enabled() is false, is_root() is true, and the
// reductions leave every density moment untouched.
void test_mpi_serial_defaults_and_noop_reduction()
{
    reset_test_config(1.0, 1.0e-6);

    expect_true("mpi_serial.enabled_is_false", !PIC::Mpi::enabled());
    expect_true("mpi_serial.is_root_is_true", PIC::Mpi::is_root());

    // Deterministic per-cell pattern so any modification is detectable.
    const auto np_of    = [](size_t i) { return 1.5 * static_cast<double>(i) + 0.25; };
    const auto up_of    = [](size_t i) { return DblVector(1.0 * i, -2.0 * i, 3.0 * i + 0.5); };
    const auto up_np_of = [](size_t i) { return DblVector(0.5 * i, 0.25 * i, -1.0 * i); };

    Grid grid(4, 4, 4, 1.0);
    for (size_t i = 0; i < grid.raw_size(); ++i) {
        Cell& cell = grid.raw_get(i);
        cell.NP    = np_of(i);
        cell.UP    = up_of(i);
        cell.UP_NP = up_np_of(i);
    }
    PIC::Mpi::allreduce_grid_density(grid);
    for (size_t i = 0; i < grid.raw_size(); ++i) {
        const Cell&       cell   = grid.raw_get(i);
        const std::string prefix = "mpi_serial.grid_cell_" + std::to_string(i);
        expect_near(prefix + ".NP", cell.NP, np_of(i), 0.0);
        expect_vector_near(prefix + ".UP", cell.UP, up_of(i), 0.0);
        expect_vector_near(prefix + ".UP_NP", cell.UP_NP, up_np_of(i), 0.0);
    }

    PIC::DensityGrid density_grid(4, 4, 4, 1.0);
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        Density& node = density_grid.raw_get(i);
        node.NP       = np_of(i);
        node.UP       = up_of(i);
        node.UP_NP    = up_np_of(i);
    }
    PIC::Mpi::allreduce_density_grid(density_grid);
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        const Density&    node   = density_grid.raw_get(i);
        const std::string prefix = "mpi_serial.density_node_" + std::to_string(i);
        expect_near(prefix + ".NP", node.NP, np_of(i), 0.0);
        expect_vector_near(prefix + ".UP", node.UP, up_of(i), 0.0);
        expect_vector_near(prefix + ".UP_NP", node.UP_NP, up_np_of(i), 0.0);
    }
}

// simulate() invokes these hooks through set_boundary_conditions(), which
// throws if a hook is missing -- the stubs make the driver runnable without a
// full case script. The group hooks receive (density_grid, group_name); the
// stubs ignore both.
void register_stub_lua_callbacks()
{
    sol::state& lua = get_lua_state();
    lua.script(R"(
        function on_iteration_begin() end
        function on_iteration_end() end
        function on_particles_moved_full_time() end
        function on_set_boundary_MF() end
        function on_set_boundary_UE() end
        function on_set_boundary_EF() end
        function on_set_boundary_eta() end
        function on_set_boundary_Te() end
        function on_set_boundary_UP() end
        function on_set_boundary_NP() end
        function on_set_boundary_group_UP(dg, name) end
        function on_set_boundary_group_NP(dg, name) end
    )");
}

// End-to-end regression for the unified driver: step-0 scatter/finalize/save,
// the 11-step main loop, the save-step group diagnostics, and the serial-mode
// MPI gates all execute. With E = B = 0 and zero resistivity the exact
// solution is inertial drift, so both the trajectory and the velocity are
// checked against closed-form values. Field finiteness guards against
// NaN/Inf escaping any stage of the merged loop.
//
// The density MUST be realistic (ne ~ 1e13 cm^-3, as in production cases):
// the Te floor (1e-4 K) creates a tiny grad(Pe) seed field at the cloud
// edge, and the Hall term UE = UP - c/(4*pi*e*ne)*curl(B) amplifies that
// seed by ~5e18/ne per step. At toy densities (ne ~ 1) the whistler CFL is
// violated and the fields explode within a few steps (observed: UE ~ 1e16,
// CFL abort at step 4); at ne ~ 1e13 the seed stays ~1e-7 G over the run
// and the drift tolerances absorb it.
void test_simulate_serial_smoke()
{
    const double h   = 1.0;
    const double tau = 1.0e-6;
    reset_test_config(h, tau);
    const int group_id = reset_single_ion_group();
    register_stub_lua_callbacks();

    const index_t n     = 16;
    const index_t steps = 5;

    auto& params                     = PIC::Config::parameters();
    params.time_steps                = steps;
    params.save_time_steps           = 1000000; // only the final step saves
    params.grid_size_x               = n;
    params.grid_size_y               = n;
    params.grid_size_z               = n;
    params.dens_cutoff               = 1.0e12;
    params.save_particle_diagnostics = false;
    params.total_particles_num       = 8;

    Grid grid(n, n, n, h);

    const DblVector v0(1.0e3, -2.0e3, 1.5e3); // cm/s; tau*v << h keeps CFL happy
    const double    ni = 1.0e13;              // realistic weight: cloud NP ~ 1e13-1e14 >> dens_cutoff

    Particles particles;
    std::vector<DblVector> starts;
    for (int dx = 0; dx < 2; ++dx)
        for (int dy = 0; dy < 2; ++dy)
            for (int dz = 0; dz < 2; ++dz) {
                const DblVector r0(7.3 + 0.4 * dx, 7.3 + 0.4 * dy, 7.3 + 0.4 * dz);
                starts.push_back(r0);
                particles.add(Particle(group_id, r0, v0, ni));
            }

    PIC::simulate(grid, particles);

    // Inertial drift: r = r0 + steps*tau*v0, v unchanged.
    const DblVector drift = static_cast<double>(steps) * tau * v0;
    for (index_t p = 0; p < particles.size(); ++p) {
        const std::string prefix   = "simulate_smoke.particle_" + std::to_string(p);
        const DblVector   expected = starts[p] + drift;
        expect_vector_near(prefix + ".r", particles[p].r, expected, 1.0e-5);
        expect_vector_near(prefix + ".v", particles[p].v, v0, 0.1, 1.0e-4);
    }
    expect_true("simulate_smoke.no_particles_lost", particles.inactive_count() == 0);

    bool all_finite = true;
    for (size_t i = 0; i < grid.raw_size(); ++i) {
        const Cell& c = grid.raw_get(i);
        all_finite = all_finite && std::isfinite(c.NP) && std::isfinite(c.Te) && std::isfinite(c.eta)
                     && std::isfinite(c.E.x) && std::isfinite(c.E.y) && std::isfinite(c.E.z)
                     && std::isfinite(c.B.x) && std::isfinite(c.B.y) && std::isfinite(c.B.z)
                     && std::isfinite(c.UE.x) && std::isfinite(c.UE.y) && std::isfinite(c.UE.z)
                     && std::isfinite(c.UP.x) && std::isfinite(c.UP.y) && std::isfinite(c.UP.z);
    }
    expect_true("simulate_smoke.grid_fields_finite", all_finite);
}

// Faraday update must keep div(B) = 0 exactly: the forward-difference
// divergence of the backward-difference curl telescopes to zero, and the
// isotropic-curl smoothing factors commute with the differences (this is the
// discrete identity the m=4-mode analysis relies on). Random E, B = 0.
void test_Faraday_div_curl_free()
{
    const double h = 1.0;
    // tau chosen so ctau_2h = c*tau/(2h) = 1: divergence residuals are then
    // pure floating-point cancellation error at O(1) field magnitude.
    reset_test_config(h, 2.0 * h / PIC::Constants::c());

    for (int iso = 0; iso <= 1; ++iso) {
        PIC::Config::parameters().isotropic_curl_enabled = (iso == 1);

        constexpr index_t n = 12;
        Grid              grid(n, n, n, h);

        std::mt19937                           gen(20260712u + iso);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (size_t i = 0; i < grid.raw_size(); ++i) {
            Cell& cell = grid.raw_get(i);
            cell.E     = DblVector(dist(gen), dist(gen), dist(gen));
            cell.B     = DblVector(0.0, 0.0, 0.0);
        }

        PIC::calc_magnetic_field_half_time(grid, PIC::Config::parameters());

        // div(B) at the integer node (i, j, k): with Bx on Edge X (i+1/2, j, k)
        // etc. the conserved divergence uses BACKWARD differences (they
        // telescope against the backward differences inside the curl), and
        // each difference D_a must carry the SAME transverse smoothing as the
        // curl: the smoothed curl conserves the smoothed divergence, not the
        // plain one (constant-coefficient stencils commute). g = 0 reduces to
        // the plain 2-point divergence for the standard curl.
        const double g  = (iso == 1) ? 1.0 / 24.0 : 0.0;
        const double w0 = 1.0 - 4.0 * g;

        const auto dx_Bx = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).B.x - grid(i - 1, j, k).B.x; };
        const auto dy_By = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).B.y - grid(i, j - 1, k).B.y; };
        const auto dz_Bz = [&grid](index_t i, index_t j, index_t k) { return grid(i, j, k).B.z - grid(i, j, k - 1).B.z; };

        // Smoothed stencil reaches +/-1 transversally and -1 along the
        // difference axis; every touched B must come from the updated
        // interior [1, n-2]: i, j, k in [2, n-3].
        double max_abs_div = 0.0;
        for (index_t i = 2; i + 2 < n; ++i)
            for (index_t j = 2; j + 2 < n; ++j)
                for (index_t k = 2; k + 2 < n; ++k) {
                    const double Dx = w0 * dx_Bx(i, j, k)
                                      + g * (dx_Bx(i, j - 1, k) + dx_Bx(i, j + 1, k)
                                             + dx_Bx(i, j, k - 1) + dx_Bx(i, j, k + 1));
                    const double Dy = w0 * dy_By(i, j, k)
                                      + g * (dy_By(i - 1, j, k) + dy_By(i + 1, j, k)
                                             + dy_By(i, j, k - 1) + dy_By(i, j, k + 1));
                    const double Dz = w0 * dz_Bz(i, j, k)
                                      + g * (dz_Bz(i - 1, j, k) + dz_Bz(i + 1, j, k)
                                             + dz_Bz(i, j - 1, k) + dz_Bz(i, j + 1, k));
                    max_abs_div = std::max(max_abs_div, std::fabs(Dx + Dy + Dz));
                }

        expect_near(std::string("Faraday_div_curl.max_abs_div.iso") + std::to_string(iso), max_abs_div, 0.0, 1.0e-12);
    }

    PIC::Config::parameters().isotropic_curl_enabled = false;
}

// CIC partition of unity: gathering a uniform field at any interior point
// must reproduce the constant exactly (weights sum to 1) for every staggering.
void test_gather_partition_of_unity()
{
    reset_test_config(1.0, 1.0e-6);

    const DblVector F0(1.7, -2.3, 0.9);
    const double    NP0 = 3.7;

    Grid grid(8, 8, 8, 1.0);
    for (size_t i = 0; i < grid.raw_size(); ++i) {
        Cell& cell = grid.raw_get(i);
        cell.E     = F0;
        cell.B     = F0;
        cell.NP    = NP0;
    }

    const std::vector<DblVector> points = {
        DblVector(3.5, 3.5, 3.5),  // cell centre
        DblVector(3.0, 3.5, 3.5),  // on a face
        DblVector(3.21, 3.47, 3.83),
        DblVector(2.5, 4.5, 3.5),
    };

    for (size_t pi = 0; pi < points.size(); ++pi) {
        const std::string prefix = "gather_unity.point_" + std::to_string(pi);

        DblVector out;
        PIC::gather_face(grid, points[pi], &Cell::E, out);
        expect_vector_near(prefix + ".face_E", out, F0, 1.0e-12, 1.0e-13);

        PIC::gather_edge(grid, points[pi], &Cell::B, out);
        expect_vector_near(prefix + ".edge_B", out, F0, 1.0e-12, 1.0e-13);

        const double np = PIC::gather_center(grid, points[pi], &Cell::NP);
        expect_near(prefix + ".center_NP", np, NP0, 1.0e-12, 1.0e-13);
    }
}

// Face-momentum identity behind set_grid_UP: the scatter deposits UP = v*w
// and UP_NP = w on the same faces, so UP/UP_NP must equal the particle
// velocity on every touched face (this is exactly the normalization
// finalize_main_grid performs).
void test_scatter_face_velocity()
{
    reset_test_config(1.0, 1.0e-12);
    const int group_id = reset_single_ion_group();

    const DblVector velocity(1.5e3, -2.5e3, 3.5e3);
    const Particle  particle(group_id, DblVector(3.37, 3.62, 3.18), velocity, 4.0);

    PIC::DensityGrid density_grid(8, 8, 8, 1.0);
    PIC::scatter_particle_snapshot(particle, density_grid);

    int touched_faces = 0;
    for (size_t i = 0; i < density_grid.raw_size(); ++i) {
        const Density& node = density_grid.raw_get(i);
        if (node.UP_NP.x > 1.0e-9) {
            expect_near("scatter_face_v.x_" + std::to_string(i), node.UP.x / node.UP_NP.x, velocity.x, 1.0e-9, 1.0e-12);
            ++touched_faces;
        }
        if (node.UP_NP.y > 1.0e-9) {
            expect_near("scatter_face_v.y_" + std::to_string(i), node.UP.y / node.UP_NP.y, velocity.y, 1.0e-9, 1.0e-12);
            ++touched_faces;
        }
        if (node.UP_NP.z > 1.0e-9) {
            expect_near("scatter_face_v.z_" + std::to_string(i), node.UP.z / node.UP_NP.z, velocity.z, 1.0e-9, 1.0e-12);
            ++touched_faces;
        }
    }
    expect_true("scatter_face_v.faces_touched", touched_faces >= 24); // 8 faces per component
}

// Thomas sweep against a dense Gaussian elimination on the same tridiagonal
// system, including the Dirichlet endpoint handling from T_old.
void test_Thomas_vs_dense()
{
    const int    K     = 12;
    const int    first = 2;
    const int    last  = K - 3;
    const double coeff = 0.35;

    std::mt19937                           gen(777);
    std::uniform_real_distribution<double> chi_dist(0.1, 3.0);
    std::uniform_real_distribution<double> te_dist(0.5, 2.0);

    std::vector<double> chi(K), T_old(K);
    for (int k = 0; k < K; ++k) {
        chi[k]   = chi_dist(gen);
        T_old[k] = te_dist(gen);
    }
    std::vector<double> T_new = T_old;

    std::vector<double> p_buf, q_buf;
    PIC::detail::run_Thomas(chi, T_old, T_new, first, last, coeff, p_buf, q_buf);

    // Dense reference: same coefficients, unknowns x[first..last].
    auto chi_face = [&](int k1, int k2) { return 2.0 * (chi[k1] * chi[k2]) / (chi[k1] + chi[k2] + 1e-20); };

    const int                        N = last - first + 1;
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double>              rhs(N, 0.0);

    for (int k = first; k <= last; ++k) {
        const int    r  = k - first;
        const double ak = -coeff * chi_face(k - 1, k);
        const double ck = -coeff * chi_face(k, k + 1);
        A[r][r]         = 1.0 - ak - ck;
        rhs[r]          = T_old[k];
        if (k > first)
            A[r][r - 1] = ak;
        else
            rhs[r] -= ak * T_old[k - 1];
        if (k < last)
            A[r][r + 1] = ck;
        else
            rhs[r] -= ck * T_old[k + 1];
    }

    // Gaussian elimination with partial pivoting.
    for (int col = 0; col < N; ++col) {
        int piv = col;
        for (int r = col + 1; r < N; ++r)
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col]))
                piv = r;
        std::swap(A[col], A[piv]);
        std::swap(rhs[col], rhs[piv]);
        for (int r = col + 1; r < N; ++r) {
            const double m = A[r][col] / A[col][col];
            for (int c2 = col; c2 < N; ++c2)
                A[r][c2] -= m * A[col][c2];
            rhs[r] -= m * rhs[col];
        }
    }
    std::vector<double> x(N);
    for (int r = N - 1; r >= 0; --r) {
        double s = rhs[r];
        for (int c2 = r + 1; c2 < N; ++c2)
            s -= A[r][c2] * x[c2];
        x[r] = s / A[r][r];
    }

    for (int r = 0; r < N; ++r)
        expect_near("Thomas_vs_dense.x_" + std::to_string(r), T_new[first + r], x[r], 1.0e-10, 1.0e-12);
}

// to_idx maps position -> home node index for each staggering offset;
// these are the exact conventions the gather/scatter stencils rely on.
void test_to_idx_centering()
{
    expect_true("to_idx.node_exact", to_idx(3.0, 1.0) == 3);
    expect_true("to_idx.below_node", to_idx(3.999999, 1.0) == 3);
    expect_true("to_idx.on_face_upper_cell", to_idx(4.0, 1.0) == 4);
    expect_true("to_idx.half_offset_below", to_idx(3.0, 1.0, 0.5) == 2);   // floor(2.5)
    expect_true("to_idx.half_offset_at_centre", to_idx(3.5, 1.0, 0.5) == 3);
    expect_true("to_idx.half_offset_just_below", to_idx(3.499999, 1.0, 0.5) == 2);
    expect_true("to_idx.fractional_h", to_idx(1.75, 0.5) == 3);
}

// Range and distribution sanity of the fixed rnd(min, max): the old formula
// u*max + min mapped rnd(2,5) onto [2,7].
void test_rnd_range()
{
    srand(20260711u);

    bool   in_bounds = true;
    double mn = 1.0e300, mx = -1.0e300, sum = 0.0;
    const int N = 10000;
    for (int i = 0; i < N; ++i) {
        const double v = rnd(2.0, 5.0);
        in_bounds      = in_bounds && v >= 2.0 && v <= 5.0;
        mn             = std::min(mn, v);
        mx             = std::max(mx, v);
        sum += v;
    }
    expect_true("rnd_range.in_bounds", in_bounds);
    expect_true("rnd_range.spans_interval", mn < 2.5 && mx > 4.5);
    expect_near("rnd_range.mean", sum / N, 3.5, 0.1);

    bool unit_in_bounds = true;
    for (int i = 0; i < 1000; ++i) {
        const double v = rnd();
        unit_in_bounds = unit_in_bounds && v >= 0.0 && v <= 1.0;
    }
    expect_true("rnd_range.default_unit", unit_in_bounds);
}

// An escaped exception (e.g. the CFL abort thrown by move_particles_half_time)
// must show up as a named FAIL, not as a silent process abort with no output.
void run_test(const char* name, void (*test_fn)())
{
    try {
        test_fn();
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << name << ": unhandled exception: " << e.what() << "\n";
        ++g_failures;
    } catch (...) {
        std::cerr << "FAIL " << name << ": unhandled non-standard exception\n";
        ++g_failures;
    }
}

} // namespace

int main()
{
    run_test("pocketfft_c2c_roundtrip", test_pocketfft_c2c_roundtrip);
    run_test("fdtd_curl_dispersion_symbol", test_fdtd_curl_dispersion_symbol);
    run_test("Boris_uniform_e", test_Boris_uniform_e);
    run_test("Boris_uniform_b_preserves_speed", test_Boris_uniform_b_preserves_speed);
    run_test("cic_scatter_weight_sum", test_cic_scatter_weight_sum);
    run_test("ngp_staggered_gather_and_scatter", test_ngp_staggered_gather_and_scatter);
    run_test("tsc_staggered_gather_and_scatter", test_tsc_staggered_gather_and_scatter);
    run_test("pstd_plane_wave_Faraday_derivative", test_pstd_plane_wave_Faraday_derivative);
    run_test("Faraday_div_curl_free", test_Faraday_div_curl_free);
    run_test("gather_partition_of_unity", test_gather_partition_of_unity);
    run_test("scatter_face_velocity", test_scatter_face_velocity);
    run_test("Thomas_vs_dense", test_Thomas_vs_dense);
    run_test("to_idx_centering", test_to_idx_centering);
    run_test("rnd_range", test_rnd_range);
    run_test("mpi_serial_defaults_and_noop_reduction", test_mpi_serial_defaults_and_noop_reduction);
    run_test("simulate_serial_smoke", test_simulate_serial_smoke);

    if (g_failures != 0) {
        std::cerr << g_failures << " OpenPIC verification test(s) failed\n";
        return 1;
    }

    std::cout << "OpenPIC verification tests passed\n";
    return 0;
}
