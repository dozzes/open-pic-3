#include "core/filters.h"

#include <vector>

namespace PIC {

// Hyperresistive induction filter: dB/dt = -nu4 * Lap(Lap(B)), where
// nu4 = eta4 * c^2 / (4*pi). Applied once per full step, AFTER the B
// corrector, and deliberately NOT as a term in Ohm's law E: the term
// exists only to drain grid-scale (2h) modes from B. Feeding it into E
// pushes particles with an unphysical force exactly where shot noise
// keeps reseeding those modes (giant cavity E spikes were observed with
// the Ohm's-law form; see docs/M20D14_Analysis_Notes.md, 2026-07-04).
// Damping ~ k^4 leaves resolved scales nearly untouched.
// Explicit stability: d4 = eta4*c^2*tau/(4*pi*h^4) <= ~0.014.
void apply_hyper_resistivity(Grid& grid, const Config::Parameters& params, SolverWorkspace& ws)
{
    const double eta4 = params.hyper_resistivity;
    if (eta4 <= 0.0)
        return;

    print_tm("apply_hyper_resistivity");

    const double inv_h2 = params.inv_h() * params.inv_h();
    const double coef   = eta4 * Constants::c() * Constants::c() / (4.0 * Constants::pi()) * params.tau;

    const index_t nx = grid.size_x();
    const index_t ny = grid.size_y();
    const index_t nz = grid.size_z();

    // Flat scratch for Lap(B); boundary entries are never written nor read
    // (first loop covers [1, n-2], the second reads only its +/-1 neighbors
    // from [2, n-3]), so no clearing is needed.
    std::vector<DblVector>& lap = ws.lap_B;
    lap.resize((size_t)nx * ny * nz);
    const auto idx = [ny, nz](index_t x, index_t y, index_t z) -> size_t { return ((size_t)x * ny + y) * nz + z; };

#pragma omp parallel for collapse(3)
    for (index_t kx = 1; kx < nx - 1; ++kx)
        for (index_t ky = 1; ky < ny - 1; ++ky)
            for (index_t kz = 1; kz < nz - 1; ++kz) {
                const DblVector& Bc = grid(kx, ky, kz).B;
                lap[idx(kx, ky, kz)] =
                    (grid(kx + 1, ky, kz).B + grid(kx - 1, ky, kz).B + grid(kx, ky + 1, kz).B
                     + grid(kx, ky - 1, kz).B + grid(kx, ky, kz + 1).B + grid(kx, ky, kz - 1).B
                     - 6.0 * Bc)
                    * inv_h2;
            }

#pragma omp parallel for collapse(3)
    for (index_t kx = 2; kx < nx - 2; ++kx)
        for (index_t ky = 2; ky < ny - 2; ++ky)
            for (index_t kz = 2; kz < nz - 2; ++kz) {
                const DblVector& Lc = lap[idx(kx, ky, kz)];
                const DblVector  lap2 =
                    (lap[idx(kx + 1, ky, kz)] + lap[idx(kx - 1, ky, kz)] + lap[idx(kx, ky + 1, kz)]
                     + lap[idx(kx, ky - 1, kz)] + lap[idx(kx, ky, kz + 1)] + lap[idx(kx, ky, kz - 1)]
                     - 6.0 * Lc)
                    * inv_h2;

                DblVector& B = grid(kx, ky, kz).B;
                B.x -= coef * lap2.x;
                B.y -= coef * lap2.y;
                B.z -= coef * lap2.z;
            }
}

} // namespace PIC
