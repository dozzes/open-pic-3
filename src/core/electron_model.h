#pragma once

#include "config/config.h"
#include "core/solver_workspace.h"
#include "grid/grid.h"

namespace PIC {

// Electron closure: resistivity model (Spitzer eta(Te) or the legacy-code-style
// kappa friction) and the electron temperature source terms. The heat
// conduction solver lives separately in core/heat_solver.h.

// Per-cell resistivity update; fills the binomially smoothed Te in the
// workspace when Spitzer_Te_smooth_passes > 0 (see smooth_te_binomial).
void update_Spitzer_coefficients(Grid& grid, const Config::Parameters& params, SolverWorkspace& workspace);

// Joule heating + adiabatic compression term for Te at cell centres.
void apply_electron_thermodynamics(Grid& grid, const Config::Parameters& params);

// Fills workspace.te_smooth with binomially smoothed Te (grid.Te untouched).
void smooth_te_binomial(const Grid& grid, SolverWorkspace& workspace, int passes);

} // namespace PIC
