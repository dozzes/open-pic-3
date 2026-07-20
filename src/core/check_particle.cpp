#include "core/check_particle.h"
#include "util/opic_fwd.h"
#include "config/config.h"
#include "io/io_utilities.h"
#include "grid/grid.h"
#include "particles/particles.h"
#include "core/gather_scatter.h"
#include "lua/call_lua_function.h"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>
#include <stdexcept>

/*********************************************************
 * Check if particle is active.                           *
 * equation of motion is solved for active particles only *
 * return true if particles is active, otherwise - false  *
 *********************************************************/
bool is_particle_can_scatter(Particles& particles, index_t p, const Grid& grid, const DblVector& dr)
{
    if (particles.is_inactive(p))
        return false;

    const Particle&      particle        = particles[p];
    const double         h               = PIC::Config::h();
    const index_t        pi              = to_idx((particle.r.x + dr.x), h);
    const index_t        pj              = to_idx((particle.r.y + dr.y), h);
    const index_t        pk              = to_idx((particle.r.z + dr.z), h);
    const PIC::CellState home_cell_state = grid(pi, pj, pk).state();

    const int kd = 1;
    if (pi < kd || pi >= (grid.size_x() - kd) || pj < kd || pj >= (grid.size_y() - kd) || pk < kd
        || pk >= (grid.size_z() - kd) || home_cell_state == PIC::cs_absorptive) {
        particles.remove_later(p);
        return false;
    }

    return true;
}

bool is_particle_can_move(Particles& particles, index_t p, const Grid& grid)
{
    if (!is_particle_can_scatter(particles, p, grid, DblVector()))
        return false;

    const double  h        = PIC::Config::h();
    Particle&     particle = particles[p];
    const index_t pi       = to_idx(particle.r.x, h);
    const index_t pj       = to_idx(particle.r.y, h);
    const index_t pk       = to_idx(particle.r.z, h);

    if (pi >= grid.size_x() || pj >= grid.size_y() || pk >= grid.size_z())
        return false;

    const PIC::CellState home_cell_state = grid(pi, pj, pk).state();
    if (home_cell_state == PIC::cs_active)
        return true;

    if (home_cell_state == PIC::cs_custom) {
        bool valid = false;

#pragma omp critical(lua_access)
        {
            valid = lua_validate_particle(particle);
        }

        return valid;
    }

    // never reach here
    return false;
}

bool check_particle_move(const Particle& particle, const Grid& grid, const DblVector& dr)
{
    if (PIC::Config::CFL_severity() == PIC::CFLSeverity::Ignore)
        return true;

    const double h   = PIC::Config::h();
    const double h_2 = PIC::Config::h_2();

    // BUG FIX: NaN/Inf in dr is NOT caught by (abs(NaN) > h_2) because that
    // comparison returns false in IEEE 754.  Must test explicitly first.
    const bool dr_nonfinite = !is_finite(dr);
    const bool dr_cfl       = std::abs(dr.x) > h_2 || std::abs(dr.y) > h_2 || std::abs(dr.z) > h_2;

    if (!dr_nonfinite && !dr_cfl)
        return true;

    // ── Diagnostics ────────────────────────────────────────────────────────
    const index_t step = PIC::Config::current_time_step();

    // BUG FIX: clamp cell indices before any grid access to avoid out-of-bounds UB.
    // particle.r may be outside the grid if it was already in bad state.
    const index_t pi = std::min(std::max(to_idx(particle.r.x, h), index_t(0)), grid.size_x() - 1);
    const index_t pj = std::min(std::max(to_idx(particle.r.y, h), index_t(0)), grid.size_y() - 1);
    const index_t pk = std::min(std::max(to_idx(particle.r.z, h), index_t(0)), grid.size_z() - 1);

    // from_grid_to_point also uses to_idx internally; only safe when r is inside domain.
    Grid::NodeType point_val;
    const bool     r_inside = particle.r.x >= 0.0 && particle.r.x < grid.size_x() * h && particle.r.y >= 0.0
                          && particle.r.y < grid.size_y() * h && particle.r.z >= 0.0 && particle.r.z < grid.size_z() * h;
    if (r_inside)
        PIC::from_grid_to_point(grid, particle.r, point_val);

    const int         tid           = omp_get_thread_num();
    const std::string log_file_name = fmt::format("opic_thread_{}_check_particle_move_err.log", tid);
    std::ofstream     ofs_log(log_file_name, std::ios_base::app);

    if (ofs_log) {
        fmt::print(ofs_log,
                   "Step = {}:\n"
                   "\t{}: particle shift dr should be < h/2 and finite:\n"
                   "\tparticle's home cell = ({}, {}, {});\n"
                   "\tgrid size = ({} x {} x {});\n"
                   "\tdr = ({}, {}, {}),\n"
                   "\tdr.abs() = {}; h/2 = {};\n"
                   "\tparticle group={} at ({}, {}, {})\n"
                   "\tparticle velocity: vx={}, vy={}, vz={}\n"
                   "\tgrid values at particle position{}:\n"
                   "\t\t B = {}, Bx = {}, By = {}, Bz = {}\n"
                   "\t\t E = {}, Ex = {}, Ey = {}, Ez = {}\n"
                   "\t\t NP = {}, UPx = {}, UPy = {}, UPz = {}\n"
                   "\t\t UEx = {}, UEy = {}, UEz = {}\n"
                   "\t\t cell state = {} (0 - active, 1 - absorptive, 2 - custom)\n",

                   step,
                   dr_nonfinite ? "NaN/Inf in dr" : "CFL violation",
                   pi,
                   pj,
                   pk,
                   grid.size_x(),
                   grid.size_y(),
                   grid.size_z(),
                   dr.x,
                   dr.y,
                   dr.z,
                   dr.abs(),
                   h_2,
                   particle.group_id,
                   particle.r.x,
                   particle.r.y,
                   particle.r.z,
                   particle.v.x,
                   particle.v.y,
                   particle.v.z, // velocity — key for diagnosis
                   r_inside ? "" : " (clamped — particle outside domain)",
                   point_val.B.abs(),
                   point_val.B.x,
                   point_val.B.y,
                   point_val.B.z,
                   point_val.E.abs(),
                   point_val.E.x,
                   point_val.E.y,
                   point_val.E.z,
                   point_val.NP,
                   point_val.UP.x,
                   point_val.UP.y,
                   point_val.UP.z,
                   point_val.UE.x,
                   point_val.UE.y,
                   point_val.UE.z,
                   static_cast<int>(grid(pi, pj, pk).state()));
    }

    return false;
}
