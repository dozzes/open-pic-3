#pragma once

#include "config/config.h"
#include "core/solver_workspace.h"
#include "grid/grid.h"
#include "io/io_utilities.h"

#include <algorithm>
#include <vector>

namespace PIC {

// Hyperresistive k^4 induction filter; see the implementation for the physics
// rationale (deliberately NOT an Ohm's-law term) and the stability bound.
void apply_hyper_resistivity(Grid& grid, const Config::Parameters& params, SolverWorkspace& workspace);

/**
 *Applies a 3D Binomial Filter(1 - 2 - 1) to the current density(UP).
 * *PHYSICS RATIONALE :
 *At high magnetolaminar parameters(delta > 10), numerical dispersion
 * on a rectangular grid causes artificial anisotropy(e.g., a spherical
 * cloud expanding as a rhombus).This filter suppresses high - frequency
 * spectral noise and restores isotropy in the equatorial plane.
 *
 *IMPLEMENTATION:
 *Separable in - place 3D filter with O(1) additional memory per thread.
 */
template <class GridT> void filter_sources(GridT& grid)
{
    print_tm("filter_sources (NP, UP, UP_NP)");

    const index_t nx = grid.size_x();
    const index_t ny = grid.size_y();
    const index_t nz = grid.size_z();

#pragma omp parallel
    {
        // Buffers for vector (UP, UP_NP) and scalar (NP) fields.
        // UP_NP MUST be smoothed with the same kernel as UP: this runs on the
        // raw accumulators BEFORE finalize_main_grid, where set_grid_UP
        // computes velocity as UP/UP_NP. Smoothing the momentum numerator but
        // not the weight denominator multiplies the velocity by
        // (smoothed w)/(raw w), which oscillates with cubic symmetry at the
        // steep cloud edge (observed as petal-shaped NP artifacts). With both
        // smoothed, UP/UP_NP is a proper weighted-average velocity.
        std::vector<DblVector> buffer_UP(std::max({nx, ny, nz}));
        std::vector<DblVector> buffer_UPNP(std::max({nx, ny, nz}));
        std::vector<double>    buffer_NP(std::max({nx, ny, nz}));

        // --- PHASE 1: Smoothing along X-axis ---
#pragma omp for
        for (index_t j = 0; j < ny; ++j) {
            for (index_t k = 0; k < nz; ++k) {
                for (index_t i = 1; i < nx - 1; ++i) {
                    buffer_UP[i]   = (grid(i - 1, j, k).UP + grid(i, j, k).UP * 2.0 + grid(i + 1, j, k).UP) * 0.25;
                    buffer_UPNP[i] = (grid(i - 1, j, k).UP_NP + grid(i, j, k).UP_NP * 2.0 + grid(i + 1, j, k).UP_NP) * 0.25;
                    buffer_NP[i]   = (grid(i - 1, j, k).NP + grid(i, j, k).NP * 2.0 + grid(i + 1, j, k).NP) * 0.25;
                }
                for (index_t i = 1; i < nx - 1; ++i) {
                    grid(i, j, k).UP    = buffer_UP[i];
                    grid(i, j, k).UP_NP = buffer_UPNP[i];
                    grid(i, j, k).NP    = buffer_NP[i];
                }
            }
        }

        // --- PHASE 2: Smoothing along Y-axis ---
#pragma omp for
        for (index_t i = 0; i < nx; ++i) {
            for (index_t k = 0; k < nz; ++k) {
                for (index_t j = 1; j < ny - 1; ++j) {
                    buffer_UP[j]   = (grid(i, j - 1, k).UP + grid(i, j, k).UP * 2.0 + grid(i, j + 1, k).UP) * 0.25;
                    buffer_UPNP[j] = (grid(i, j - 1, k).UP_NP + grid(i, j, k).UP_NP * 2.0 + grid(i, j + 1, k).UP_NP) * 0.25;
                    buffer_NP[j]   = (grid(i, j - 1, k).NP + grid(i, j, k).NP * 2.0 + grid(i, j + 1, k).NP) * 0.25;
                }
                for (index_t j = 1; j < ny - 1; ++j) {
                    grid(i, j, k).UP    = buffer_UP[j];
                    grid(i, j, k).UP_NP = buffer_UPNP[j];
                    grid(i, j, k).NP    = buffer_NP[j];
                }
            }
        }

        // --- PHASE 3: Smoothing along Z-axis ---
#pragma omp for
        for (index_t i = 0; i < nx; ++i) {
            for (index_t j = 0; j < ny; ++j) {
                for (index_t k = 1; k < nz - 1; ++k) {
                    buffer_UP[k]   = (grid(i, j, k - 1).UP + grid(i, j, k).UP * 2.0 + grid(i, j, k + 1).UP) * 0.25;
                    buffer_UPNP[k] = (grid(i, j, k - 1).UP_NP + grid(i, j, k).UP_NP * 2.0 + grid(i, j, k + 1).UP_NP) * 0.25;
                    buffer_NP[k]   = (grid(i, j, k - 1).NP + grid(i, j, k).NP * 2.0 + grid(i, j, k + 1).NP) * 0.25;
                }
                for (index_t k = 1; k < nz - 1; ++k) {
                    grid(i, j, k).UP    = buffer_UP[k];
                    grid(i, j, k).UP_NP = buffer_UPNP[k];
                    grid(i, j, k).NP    = buffer_NP[k];
                }
            }
        }
    }
}

} // namespace PIC
