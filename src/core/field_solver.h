#pragma once

#include "config/config.h"
#include "core/solver_workspace.h"
#include "grid/grid.h"

namespace PIC {

// Electromagnetic field solver on the Yee (C-grid) stencil.
// All functions take the run parameters explicitly: they depend on global
// state only through the caller-supplied Parameters.

// Faraday half-step B(n) -> B(n+1/2), finite-difference curl (production).
void calc_magnetic_field_half_time(Grid& grid, const Config::Parameters& params);

// Faraday half-step with a spectral (FFT) curl. Experimental: non-periodic
// boundaries make it unusable for production physics (see PSTD status in the
// User Guide); kept for method development.
void calc_magnetic_field_pstd_half_time(Grid& grid, const Config::Parameters& params, SolverWorkspace& workspace);

// Electron velocity from Ampere's law: Ue = Ui - curl(B)/(4*pi*e*n).
void calc_electrons_velocity(Grid& grid, const Config::Parameters& params);

// Generalised Ohm's law: E = (UE x B)/c - grad(Pe)/(e*ne) + eta*(c/4pi)*curl(B).
void calc_electric_field(Grid& grid, const Config::Parameters& params);

} // namespace PIC
