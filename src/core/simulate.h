#pragma once

// Per-step functions for the hybrid PIC advance and the top-level simulation loop.
//
// Time-step cycle (leapfrog):
//   I.   B(m) -> B(m+1/2)           calc_magnetic_field_half_time (core/field_solver.h)
//   II.  push particles to (m+1/2)  move_particles_half_time  +  scatter -> grid
//   III. scatter -> grid             (see gather_scatter.h: from_particles_to_grid)
//   IV.  push particles to (m+1)    move_particles_full_time
//   V.   B(m+1/2) -> B(m+1)         calc_magnetic_field_half_time (repeated)
//   VI.  UE, E from Ohm's law        calc_electrons_velocity, calc_electric_field
//
// The field solver lives in core/field_solver.h, the electron closure in
// core/electron_model.h, grid filters in core/filters.h; this header keeps
// the particle-advance functions and the driver.

#include "util/opic_fwd.h"
#include "grid/grid.h"
#include "core/field_solver.h"
#include "core/solver_workspace.h"
#include <string>

namespace PIC {

typedef GridContainer<Density> DensityGrid;
void                           move_particles_half_time(std::vector<DensityGrid>& density_grids,
                                                        const Grid&               grid,
                                                        Particles&                particles,
                                                        const std::string&        group_name);

void move_particles_full_time(const Grid& grid, Particles& particles, const std::string& group_name);

void push_particle_Boris(const Grid& grid, Particle& particle);
void scatter_particle_snapshot(const Particle& particle, DensityGrid& grid);

// Single driver for both serial and MPI replicated-grid runs; the mode is
// taken from Mpi::enabled(). In MPI mode each rank holds the full grid and a
// slice of the particles: after each scatter step all ranks sum their density
// moments via MPI_Allreduce (a no-op in serial runs), file output happens on
// the root rank only, and per-particle diagnostics are skipped (particles are
// rank-distributed and there is no gather step).
void simulate(Grid& grid, Particles& particles);

} // namespace PIC
