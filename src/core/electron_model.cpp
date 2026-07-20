#include "core/electron_model.h"
#include "core/cell_stencil.h"
#include "io/io_utilities.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace PIC {

// Fills `ws.te_smooth` with Te smoothed by `passes` separable 3-tap binomial
// ([1,2,1]/4) passes (x, then y, then z per pass). A single pass exactly
// cancels a period-2h checkerboard: at a checkerboard node v(i-1)=v(i+1)=-v(i),
// so (v(i-1)+2*v(i)+v(i+1))/4 = 0. Structure spanning several cells (e.g. the
// compressed front, 3+ cells wide) passes through with little attenuation.
// Boundary cells (no interior neighbor along that axis) keep their value.
// `grid.Te` itself is never modified -- this only prepares the value fed to
// the eta(Te) formula below.
void smooth_te_binomial(const Grid& grid, SolverWorkspace& ws, int passes)
{
    const index_t nx = grid.size_x();
    const index_t ny = grid.size_y();
    const index_t nz = grid.size_z();
    const size_t  n  = (size_t)nx * (size_t)ny * (size_t)nz;

    auto idx = [ny, nz](index_t x, index_t y, index_t z) -> size_t {
        return ((size_t)x * ny + y) * nz + z;
    };

    std::vector<double>& a = ws.te_pass_a;
    std::vector<double>& b = ws.te_pass_b;
    a.resize(n);
    b.resize(n);

#pragma omp parallel for
    for (index_t x = 0; x < nx; ++x)
        for (index_t y = 0; y < ny; ++y)
            for (index_t z = 0; z < nz; ++z)
                a[idx(x, y, z)] = grid(x, y, z).Te;

    for (int pass = 0; pass < passes; ++pass) {
#pragma omp parallel for
        for (index_t x = 0; x < nx; ++x)
            for (index_t y = 0; y < ny; ++y)
                for (index_t z = 0; z < nz; ++z) {
                    const size_t i = idx(x, y, z);
                    b[i]           = (x > 0 && x < nx - 1)
                                        ? 0.25 * (a[idx(x - 1, y, z)] + 2.0 * a[i] + a[idx(x + 1, y, z)])
                                        : a[i];
                }

#pragma omp parallel for
        for (index_t x = 0; x < nx; ++x)
            for (index_t y = 0; y < ny; ++y)
                for (index_t z = 0; z < nz; ++z) {
                    const size_t i = idx(x, y, z);
                    a[i]           = (y > 0 && y < ny - 1)
                                        ? 0.25 * (b[idx(x, y - 1, z)] + 2.0 * b[i] + b[idx(x, y + 1, z)])
                                        : b[i];
                }

#pragma omp parallel for
        for (index_t x = 0; x < nx; ++x)
            for (index_t y = 0; y < ny; ++y)
                for (index_t z = 0; z < nz; ++z) {
                    const size_t i = idx(x, y, z);
                    b[i]           = (z > 0 && z < nz - 1)
                                        ? 0.25 * (a[idx(x, y, z - 1)] + 2.0 * a[i] + a[idx(x, y, z + 1)])
                                        : a[i];
                }

        a.swap(b);
    }

    ws.te_smooth.swap(a);
}

// Spitzer resistivity: eta(Te) = eta_ref * (Te_ref / Te)^1.5
//
// eta_ref (= Config::resistivity) is the physical Spitzer value at Te_ref [s].
// The ratio form is dimensionally consistent: eta [s], Te_ref/Te [dimensionless].
// When Te_ref == 0 the scaling is disabled and eta = eta_ref everywhere.
//
// Cells below dens_cutoff have been filled artificially by set_threshold.
// Te is unphysical there and must be pinned to Te_ref to prevent runaway
// via dTe ~ eta * j^2 / NP with NP at its artificial minimum.
//
// The same runaway can start one cell OUTSIDE the floored region: a cell
// that sits just above dens_cutoff (so it is NOT pinned) but is face-
// adjacent to a floored neighbor inherits that neighbor's artificial Te/UE
// discontinuity through curl(B) and grad(Pe) -- both computed from that
// neighbor's field values -- without inheriting the neighbor's protective
// skip. With eta pinned at eta_floor (FRONTFLOOR: eta_floor=base_eta, never
// collapsing toward zero) and no conduction sink (chi=0 in some cases), the
// resulting Joule heating at this interface cell has no bound and no way to
// diffuse away; observed in practice as a step-by-step Te runaway in a
// single background-lattice corner cell (docs/WORK_LOG.md, 2026-07-18 entry:
// MA20_D14_UNIF_TE_FRONTFLOOR_HFINE_BIGBOX crash at step 363). Both the
// Te-pin below and the heating-skip in apply_electron_thermodynamics
// therefore also cover any cell with a floored face-neighbor, not only
// floored cells themselves.
//
// Spitzer_Te_smooth_passes > 0 feeds a binomial-smoothed Te into the eta
// ratio instead of the raw cell.Te (see smooth_te_binomial above). This
// targets the specific noise -> eta(Te) -> grid-mode feedback loop
// (docs/M20D14_Analysis_Notes.md) without disabling genuine large-scale
// Te-dependence the way Spitzer_Te_ref=0 (ETACONST) does. Measured (Ma=20.3,
// 2026-07-05 entry) to have no effect: the eta(r) trend across the
// cavity/front is smooth over tens of cells, not grid-scale noise, so a
// single-cell filter doesn't touch it.
//
// Spitzer_eta_floor_mult raises eta_floor from its permissive default
// (base_eta*1e-3) toward base_eta itself. The same measurement showed the
// adaptive formula's actual behavior is two-sided: eta pins at eta_ceil
// (~2.5x base_eta) in the cold cavity where Te < Te_ref, but drops to
// ~0.5-0.8x base_eta in the hot compressed shell/front where Te > Te_ref --
// exactly where the m=4 mode grows. Setting Spitzer_eta_floor_mult=1.0
// raises eta at the front back up to base_eta (like ETACONST there) while
// leaving the cold cavity's physically-adaptive ceiling behavior untouched.
//
// kappa_friction > 0 (with kappa_B_ref > 0) REPLACES the Spitzer model with
// the legacy code's constant-friction closure: fielde adds xappa*(Ui-Ue) to E
// with xappa a fixed dimensionless scalar, and its Ue = Ui - curl(H)/n
// (block3), so the resistive term is xappa*curl(B*)/n* -- divided by local
// density. Denormalized (E* = E/e0, U* = U/va, e0/va = B0/c, Ui-Ue = j/(e*ne))
// this is a per-cell resistivity
//     eta = kappa_friction * kappa_B_ref / (c * e * ne)   [s]
// with NO Te feedback and a 1/ne profile: weak in the dense shell/core,
// strong in rarefied regions. Only the explicit-diffusion stability cap
// (eta_stab) applies -- the Spitzer floor/ceiling are base_eta multiples
// and don't belong to this mode. The same cell.eta then feeds the Joule
// term of apply_electron_thermodynamics, reproducing the legacy code's temp2
// friction heating 0.5*tau*xappa*|Ui-Ue|^2.
void update_Spitzer_coefficients(Grid& grid, const Config::Parameters& params, SolverWorkspace& ws)
{
    print_tm("update_Spitzer_coefficients");

    const double base_eta      = params.resistivity;
    const double Te_ref        = params.Spitzer_Te_ref;
    const int    smooth_passes = params.Spitzer_Te_smooth_passes;
    const double eta_floor     = base_eta * params.Spitzer_eta_floor_mult;
    // The resistive term in Ohm's law is integrated explicitly, so the
    // diffusion number D = eta*c^2*tau/(4*pi*h^2) must stay below ~1/6.
    // In adiabatically cooled cloud cores eta ~ Te^-1.5 otherwise grows
    // past that limit and drives a grid-scale checkerboard in B (see
    // docs/M20D14_Analysis_Notes.md, 2026-07-02/03 entries). The ceiling
    // is therefore tied to the stability criterion itself (D <= 0.1),
    // not to a fixed multiple of base_eta: a fixed 2.5x cap chosen at
    // Ma=20 parameters (D0 ~ 0.04) already over-admits at M2D43, where
    // tau is twice as large (D0 ~ 0.083).
    const double h_grid    = params.h;
    const double tau_step  = params.tau;
    const double eta_stab  = 0.1 * 4.0 * Constants::pi() * h_grid * h_grid
                             / (Constants::c() * Constants::c() * tau_step);
    const double eta_ceil  = std::min(base_eta * 2.5, eta_stab);

    // Legacy-code xappa mode (see the function comment): eta = kappa*B0/(c*e*ne).
    const double kappa       = params.kappa_friction;
    const double kappa_B     = params.kappa_B_ref;
    const bool   use_kappa   = kappa > 0.0 && kappa_B > 0.0;
    const double kappa_coeff = use_kappa
                                   ? kappa * kappa_B / (Constants::c() * Constants::e())
                                   : 0.0;

    const index_t nx = grid.size_x();
    const index_t ny = grid.size_y();
    const index_t nz = grid.size_z();

    const bool use_smooth = Te_ref > 0.0 && smooth_passes > 0;
    if (use_smooth)
        smooth_te_binomial(grid, ws, smooth_passes);

    // True if this cell or any face-adjacent neighbor is below dens_cutoff
    // (see the function comment above). Cells at the grid boundary simply
    // have fewer neighbors to check.
    const double dens_cutoff = params.dens_cutoff;
    auto         near_floor  = [&](index_t x, index_t y, index_t z) -> bool {
        if (grid(x, y, z).NP < dens_cutoff)
            return true;
        if (x > 0 && grid(x - 1, y, z).NP < dens_cutoff)
            return true;
        if (x < nx - 1 && grid(x + 1, y, z).NP < dens_cutoff)
            return true;
        if (y > 0 && grid(x, y - 1, z).NP < dens_cutoff)
            return true;
        if (y < ny - 1 && grid(x, y + 1, z).NP < dens_cutoff)
            return true;
        if (z > 0 && grid(x, y, z - 1).NP < dens_cutoff)
            return true;
        if (z < nz - 1 && grid(x, y, z + 1).NP < dens_cutoff)
            return true;
        return false;
    };

#pragma omp parallel for collapse(3)
    for (index_t x = 0; x < nx; ++x)
        for (index_t y = 0; y < ny; ++y)
            for (index_t z = 0; z < nz; ++z) {
                Cell&        cell   = grid(x, y, z);
                const size_t i      = ((size_t)x * ny + y) * nz + z;
                const double eta_te = use_smooth ? ws.te_smooth[i] : cell.Te;

                if (use_kappa) {
                    // Constant dimensionless friction, no Te feedback; 1/ne diverges
                    // toward the cavity/vacuum, so the explicit-diffusion stability
                    // cap is the one clamp this mode needs.
                    const double ne = std::max(cell.NP, params.dens_cutoff);
                    cell.eta        = std::min(kappa_coeff / ne, eta_stab);
                    continue;
                }

                if (Te_ref > 0.0 && cell.Te > 1e-2)
                    cell.eta = base_eta * std::pow(Te_ref / eta_te, 1.5);
                else
                    cell.eta = base_eta;

                // Pin vacuum cells (and cells adjacent to a vacuum cell) to
                // reference state.
                if (Te_ref > 0.0 && near_floor(x, y, z)) {
                    cell.Te  = Te_ref;
                    cell.eta = base_eta;
                }

                if (cell.eta < eta_floor)
                    cell.eta = eta_floor;
                if (cell.eta > eta_ceil)
                    cell.eta = eta_ceil;
            }
}

// Joule: dTe = eta*|j|²*tau / (1.5*n*k_B)
// Adiabatic: dTe = -(2/3)*Te*div(Ue)*tau
// Both terms at cell centre (i+½,j+½,k+½).
// j = (c/4π)*curl(B) lives on faces; averaging the two bracketing faces
// brings it to the cell centre without widening the stencil.
// div(Ue): one-cell forward difference — exact for C-grid topology.
void apply_electron_thermodynamics(Grid& grid, const Config::Parameters& params)
{
    constexpr double c_4pi           = Constants::c_4pi();
    constexpr double k_B             = Constants::k_B();
    const double     inv_h           = 1.0 / grid.step();
    const double     joule_coeff     = params.tau / (1.5 * k_B);
    const double     adiabatic_coeff = (2.0 / 3.0) * params.tau;
    const double     dens_cutoff     = params.dens_cutoff;

    const index_t to_kx = grid.size_x() - 1;
    const index_t to_ky = grid.size_y() - 1;
    const index_t to_kz = grid.size_z() - 1;

#pragma omp parallel for collapse(3)
    for (index_t kx = 1; kx < to_kx; ++kx)
        for (index_t ky = 1; ky < to_ky; ++ky)
            for (index_t kz = 1; kz < to_kz; ++kz) {
                const CellStencil s {grid, kx, ky, kz};
                Cell&             cell = s.c();

                // Cells below dens_cutoff were filled artificially by set_threshold;
                // skip them, and any cell face-adjacent to one, to avoid runaway
                // Joule heating dTe ~ eta*j^2/NP at the floor/real-density interface
                // (see the function comment above update_Spitzer_coefficients).
                if (cell.NP < dens_cutoff || s.xm().NP < dens_cutoff || s.xp().NP < dens_cutoff
                    || s.ym().NP < dens_cutoff || s.yp().NP < dens_cutoff || s.zm().NP < dens_cutoff
                    || s.zp().NP < dens_cutoff)
                    continue;

                // ── 1. Joule heating: dTe = eta * |j|^2 * tau / (1.5 * n * k_B) ──────
                //
                // j = (c/4π) * curl(B) at cell centre (i+½, j+½, k+½).
                // curl(B)_x lives on Face X; average the two bracketing faces to
                // reach the cell centre.  Same logic applies to jy and jz.
                const double half_c4pi_invh = 0.5 * c_4pi * inv_h;

                const double jx = half_c4pi_invh
                                  * (+(s.yp().B.z - cell.B.z) + (s.xp_yp().B.z - s.xp().B.z)
                                     - (s.zp().B.y - cell.B.y) - (s.xp_zp().B.y - s.xp().B.y));

                const double jy = half_c4pi_invh
                                  * (+(s.zp().B.x - cell.B.x) + (s.yp_zp().B.x - s.yp().B.x)
                                     - (s.xp().B.z - cell.B.z) - (s.xp_yp().B.z - s.yp().B.z));

                const double jz = half_c4pi_invh
                                  * (+(s.xp().B.y - cell.B.y) + (s.xp_zp().B.y - s.zp().B.y)
                                     - (s.yp().B.x - cell.B.x) - (s.yp_zp().B.x - s.zp().B.x));

                const double dTe_joule = cell.eta * (jx * jx + jy * jy + jz * jz) * joule_coeff / cell.NP;

                // ── 2. Adiabatic cooling/heating: dTe = -(2/3) * Te * div(Ue) * tau ──
                //
                // UEx lives on Face X at integer x; forward difference over one cell
                // is the exact C-grid divergence centred at (i+½, j+½, k+½).
                const double div_Ue = inv_h
                                      * ((s.xp().UE.x - cell.UE.x) + (s.yp().UE.y - cell.UE.y)
                                         + (s.zp().UE.z - cell.UE.z));

                const double dTe_adiabatic = -cell.Te * div_Ue * adiabatic_coeff;

                // ── 3. Update ────────────────────────────────────────────────────────
                cell.Te = std::max(cell.Te + dTe_joule + dTe_adiabatic, 1e-4);
            }
}

} // namespace PIC
