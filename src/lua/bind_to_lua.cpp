#include "lua/bind_to_lua.h"
#include "config/config.h"
#include "lua/use_lua.h"
#include "core/gather_scatter.h"
#include "grid/grid.h"
#include "grid/grid_filters.h"
#include "particles/particles.h"
#include "particles/particle_groups.h"
#include "particles/marker_particles.h"
#include "io/save_grid.h"

#include <fmt/core.h>
#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>
#include <string>
#include <iostream>
#include <stdexcept>

namespace {

// --- LUA ERROR HANDLING ---
void check_sol_result(const sol::protected_function_result& result, const char* filename)
{
    if (!result.valid()) {
        sol::error err = result;
        std::cerr << "\nLua Error in \"" << filename << "\":\n" << err.what() << std::endl;
        throw std::runtime_error("Lua script execution failed.");
    }
}

} // unnamed namespace

// --- THERMAL INITIALIZATION ---
void init_thermal_from_lua(Grid& grid, sol::state& lua)
{
    auto temp_func = lua["get_initial_Te"];

    for (index_t kx = 0; kx < grid.size_x(); ++kx) {
        for (index_t ky = 0; ky < grid.size_y(); ++ky) {
            for (index_t kz = 0; kz < grid.size_z(); ++kz) {
                // Get coordinates for Lua logic
                double x = kx * grid.step();
                double y = ky * grid.step();
                double z = kz * grid.step();

                // Call Lua function for each cell
                grid(kx, ky, kz).Te = temp_func(x, y, z);
            }
        }
    }
}

bool bind_to_lua(sol::state& lua, const char* lua_cfg_file_name, Grid& grid, Particles& particles)
{
    using namespace PIC;

    // --- INITIALIZE LUA LIBRARIES ---
    lua.open_libraries(sol::lib::base,
                       sol::lib::os,
                       sol::lib::io,
                       sol::lib::table,
                       sol::lib::package,
                       sol::lib::debug,
                       sol::lib::string,
                       sol::lib::math);
    lua["package"]["path"] = lua["package"]["path"].get<std::string>() + ";./?.lua";

    // Lua 5.4 removed math.pow; task scripts and sim/lib predate that.
    lua.script("math.pow = math.pow or function(x, y) return x ^ y end");

    // --- BIND GLOBAL FUNCTIONS ---
    lua.set_function("pic_gather", sol::overload([](const Grid& g, const DblVector& v, Grid::NodeType& n) {
                         from_grid_to_point(g, v, n);
                     }));
    lua.set_function("pic_save_grid_node", &save_grid_node);

    // --- BIND DblVector ---
    lua.new_usertype<DblVector>(
        "DblVector",
        sol::call_constructor,
        sol::constructors<DblVector(), DblVector(const DblVector&), DblVector(double, double, double)>(),
        "x",
        sol::property(&DblVector::get_x, &DblVector::set_x),
        "y",
        sol::property(&DblVector::get_y, &DblVector::set_y),
        "z",
        sol::property(&DblVector::get_z, &DblVector::set_z),
        "abs",
        &DblVector::abs);

    // --- BIND ENUMS ---
    lua.new_enum("CellState",
                 "cs_active",
                 CellState::cs_active,
                 "cs_absorptive",
                 CellState::cs_absorptive,
                 "cs_custom",
                 CellState::cs_custom);

    lua.new_enum("Diagnostics",
                 "no_diag",
                 Diagnostics::no_diag,
                 "save_positions",
                 Diagnostics::save_positions,
                 "save_grid_values",
                 Diagnostics::save_grid_values,
                 "save_all",
                 Diagnostics::save_all);

    // --- BIND CELL ---
    lua.new_usertype<Cell>(
        "Cell",
        sol::constructors<Cell(),
                          Cell(double, const DblVector&, const DblVector&, const DblVector&, const DblVector&, double, double&),
                          Cell(const Cell&)>(),
        "NP",
        &Cell::NP,
        "B",
        &Cell::B,
        "E",
        &Cell::E,
        "UE",
        &Cell::UE,
        "UP",
        &Cell::UP,
        "Te",
        &Cell::Te,
        "eta",
        &Cell::eta,
        "state",
        sol::property(&Cell::state, &Cell::set_state));

    // --- BIND GRID ---
    lua.new_usertype<Grid>("Grid",
                           "at",
                           sol::overload(static_cast<Cell& (Grid::*)(index_t, index_t, index_t)>(&Grid::at),
                                         static_cast<const Cell& (Grid::*)(index_t, index_t, index_t) const>(&Grid::at)),
                           "set",
                           &Grid::set,
                           "size_x",
                           &Grid::size_x,
                           "size_y",
                           &Grid::size_y,
                           "size_z",
                           &Grid::size_z,
                           "resize",
                           &Grid::resize,
                           "length_x",
                           &Grid::length_x,
                           "length_y",
                           &Grid::length_y,
                           "length_z",
                           &Grid::length_z,
                           "set_boundary_state",
                           &Grid::set_boundary_state,
                           "step",
                           sol::property(&Grid::step, &Grid::set_step));

    // --- BIND DENSITY ---
    lua.new_usertype<Density>(
        "Density", sol::call_constructor, sol::constructors<Density()>(), "UP", &Density::UP, "NP", &Density::NP);

    // --- BIND DENSITYGRID ---
    lua.new_usertype<DensityGrid>(
        "DensityGrid",
        "at",
        sol::overload(static_cast<Density& (DensityGrid::*)(index_t, index_t, index_t)>(&DensityGrid::at),
                      static_cast<const Density& (DensityGrid::*)(index_t, index_t, index_t) const>(&DensityGrid::at)),
        "size_x",
        &DensityGrid::size_x,
        "size_y",
        &DensityGrid::size_y,
        "size_z",
        &DensityGrid::size_z,
        "step",
        sol::property(&DensityGrid::step, &DensityGrid::set_step),
        "resize",
        &DensityGrid::resize,
        "cell_volume",
        &DensityGrid::cell_volume);

    // --- BIND USERGRIDFILTERS ---
    lua.new_usertype<UserGridFilters>("GridSaveFilters", "at", &UserGridFilters::at, "add", &UserGridFilters::append);

    // --- BIND PARTICLE ---
    lua.new_usertype<Particle>("Particle",
                               sol::constructors<Particle(),
                                                 Particle(const Particle&),
                                                 Particle(int, const DblVector&, const DblVector&, double),
                                                 Particle(int, double, double, double, double, double, double, double)>(),
                               "ni",
                               &Particle::ni,
                               "r",
                               &Particle::r,
                               "v",
                               &Particle::v,
                               "group_id",
                               sol::readonly(&Particle::group_id));

    // --- BIND PARTICLES ---
    lua.new_usertype<Particles>(
        "Particles",
        sol::constructors<Particles()>(),
        "at",
        sol::overload(static_cast<Particle& (Particles::*)(index_t)>(&Particles::at),
                      static_cast<const Particle& (Particles::*)(index_t) const>(&Particles::at)),
        "set",
        sol::overload(
            static_cast<void (Particles::*)(index_t, const Particle&)>(&Particles::set),
            static_cast<void (Particles::*)(index_t, const std::string&, const DblVector&, const DblVector&, double)>(
                &Particles::set),
            static_cast<void (Particles::*)(index_t, int, double, double, double, double, double, double, double)>(
                &Particles::set)),
        "add",
        sol::overload(static_cast<void (Particles::*)(const Particle&)>(&Particles::add),
                      static_cast<void (Particles::*)(const std::string&, const DblVector&, const DblVector&, double)>(
                          &Particles::add),
                      static_cast<void (Particles::*)(int, double, double, double, double, double, double, double)>(
                          &Particles::add)),
        "erase",
        &Particles::erase,
        "size",
        sol::property(&Particles::size, &Particles::resize));

    // --- BIND PARTICLEGROUPS ---
    lua.new_usertype<ParticleGroups::ParticleGroup>("ParticleGroup",
                                                    "name",
                                                    &ParticleGroups::ParticleGroup::name,
                                                    "id",
                                                    &ParticleGroups::ParticleGroup::id,
                                                    "charge",
                                                    &ParticleGroups::ParticleGroup::charge,
                                                    "mass",
                                                    &ParticleGroups::ParticleGroup::mass,
                                                    "boundary_kind",
                                                    &ParticleGroups::ParticleGroup::boundary_kind);

    lua.new_usertype<ParticleGroups>("ParticleGroups",
                                     "create_group",
                                     &ParticleGroups::create_group,
                                     "set_boundary_kind",
                                     &ParticleGroups::set_boundary_kind,
                                     "get_group",
                                     sol::overload([](ParticleGroups&    self,
                                                      const std::string& name) { return self.get_group(name); },
                                                   [](ParticleGroups& self, int id) { return self.get_group(id); }));

    // --- BIND MARKERPARTICLES ---
    lua.new_usertype<MarkerParticles>("MarkerParticles", "add", &MarkerParticles::insert_idx);

    // --- BIND CONSTANTS ---
    sol::table pic       = lua.create_named_table("PIC");
    sol::table constants = pic.create_named("Constants");

    constants["c"]  = PIC::Constants::c();
    constants["e"]  = PIC::Constants::e();
    constants["mp"] = PIC::Constants::mp();
    constants["me"] = PIC::Constants::me();

    // --- BIND PARAMETER ENUMS ---
    lua.new_enum("ParticlePushAlg", "Direct", ParticlePushAlg::Direct, "Boris", ParticlePushAlg::Boris);

    lua.new_enum("ScatterAlg", "Standard", ScatterAlg::Standard, "NGP", ScatterAlg::NGP, "TSC", ScatterAlg::TSC
                 // "Zigzag",    ScatterAlg::Zigzag,    // disabled: no UP_NP accumulation
                 // "Esirkepov", ScatterAlg::Esirkepov  // disabled: no UP_NP accumulation
    );

    lua.new_enum("GridThreshold", "Min_Density", GridThreshold::Min_Density, "Local_CFL", GridThreshold::Local_CFL);

    lua.new_enum("CFLSeverity", "Ignore", CFLSeverity::Ignore, "Absorb", CFLSeverity::Absorb, "Stop", CFLSeverity::Stop);

    lua.new_enum("MagneticFieldAlg", "FDTD", MagneticFieldAlg::FDTD, "PSTD", MagneticFieldAlg::PSTD);

    lua.new_enum("BoundaryKind",
                 "Default",
                 ParticleGroups::BoundaryKind::Default,
                 "FlowBackground",
                 ParticleGroups::BoundaryKind::FlowBackground);

    // --- BIND CONFIG PARAMETERS ---
    lua.new_usertype<Config::Parameters>("Parameters",
                                         "tau",
                                         &Config::Parameters::tau,
                                         "time_steps",
                                         &Config::Parameters::time_steps,
                                         "current_time_step",
                                         sol::readonly(&Config::Parameters::current_time_step),
                                         "save_time_steps",
                                         &Config::Parameters::save_time_steps,
                                         "dens_cutoff",
                                         &Config::Parameters::dens_cutoff,
                                         "save_all_particles",
                                         &Config::Parameters::save_all_particles,
                                         "save_particle_diagnostics",
                                         &Config::Parameters::save_particle_diagnostics,
                                         "save_raw_grid_debug",
                                         &Config::Parameters::save_raw_grid_debug,
                                         "save_whole_grid",
                                         &Config::Parameters::save_whole_grid,
                                         "save_grid_x_plains",
                                         &Config::Parameters::save_grid_x_plains,
                                         "save_grid_y_plains",
                                         &Config::Parameters::save_grid_y_plains,
                                         "save_grid_z_plains",
                                         &Config::Parameters::save_grid_z_plains,
                                         "os_name",
                                         sol::readonly(&Config::Parameters::os_name),
                                         "process_idx",
                                         sol::readonly(&Config::Parameters::process_rank),
                                         "process_num",
                                         sol::readonly(&Config::Parameters::process_num),
                                         "L_scale",
                                         &Config::Parameters::L_scale,
                                         "T_scale",
                                         &Config::Parameters::T_scale,
                                         "U_scale",
                                         &Config::Parameters::U_scale,
                                         "N_scale",
                                         &Config::Parameters::N_scale,
                                         "E_scale",
                                         &Config::Parameters::E_scale,
                                         "B_scale",
                                         &Config::Parameters::B_scale,
                                         "push_method",
                                         &Config::Parameters::particle_push_alg,
                                         "scatter_method",
                                         &Config::Parameters::scatter_alg,
                                         "grid_threshold",
                                         &Config::Parameters::grid_threshold,
                                         "electron_thermal_conductivity",
                                         &Config::Parameters::electron_thermal_conductivity,
                                         "density_threshold",
                                         &Config::Parameters::grid_threshold,
                                         "cfl_severity",
                                         &Config::Parameters::CFL_severity,
                                         "resistivity",
                                         &Config::Parameters::resistivity,
                                         "Spitzer_Te_ref",
                                         &Config::Parameters::Spitzer_Te_ref,
                                         "Spitzer_Te_smooth_passes",
                                         &Config::Parameters::Spitzer_Te_smooth_passes,
                                         "Spitzer_eta_floor_mult",
                                         &Config::Parameters::Spitzer_eta_floor_mult,
                                         "hyper_resistivity",
                                         &Config::Parameters::hyper_resistivity,
                                         "kappa_friction",
                                         &Config::Parameters::kappa_friction,
                                         "kappa_B_ref",
                                         &Config::Parameters::kappa_B_ref,
                                         "use_filtering",
                                         &Config::Parameters::use_filtering,
                                         "random_seed",
                                         &Config::Parameters::random_seed,
                                         "magnetic_field_alg",
                                         &Config::Parameters::magnetic_field_alg,
                                         "isotropic_curl_enabled",
                                         &Config::Parameters::isotropic_curl_enabled,
                                         "cold_electrons_enabled",
                                         &Config::Parameters::cold_electrons_enabled);

    try {
        // --- ASSIGN GLOBAL INSTANCES ---
        lua["pic_grid"]      = &grid;
        lua["pic_particles"] = &particles;

        static UserGridFilters grid_save_filters;
        lua["pic_grid_save_filters"] = &grid_save_filters;

        static ParticleGroups part_groups;
        lua["pic_particle_groups"] = &part_groups;

        static MarkerParticles markers;
        lua["pic_marker_particles"] = &markers;

        Config::Parameters& parameters = Config::parameters();
        parameters.cfg_script_name     = lua_cfg_file_name;
        lua["pic_parameters"]          = &parameters;

        // --- CONFIG SETTER METHODS (called from main.lua before sim_core) ---
        // Usage: pic_config.set_Spitzer_eta_floor_mult(1.0) etc.
        sol::table pic_config = lua.create_table();
        pic_config.set_function("set_Spitzer_eta_floor_mult",
            [&parameters](double val) { parameters.set_Spitzer_eta_floor_mult(val); });
        pic_config.set_function("set_Spitzer_Te_ref",
            [&parameters](double val) { parameters.set_Spitzer_Te_ref(val); });
        pic_config.set_function("set_isotropic_curl_enabled",
            [&parameters](bool val) { parameters.set_isotropic_curl_enabled(val); });
        pic_config.set_function("set_cold_electrons_enabled",
            [&parameters](bool val) { parameters.set_cold_electrons_enabled(val); });
        pic_config.set_function("set_kappa_friction",
            [&parameters](double val) { parameters.set_kappa_friction(val); });
        pic_config.set_function("set_resistivity",
            [&parameters](double val) { parameters.set_resistivity(val); });
        pic_config.set_function("set_verbosity_level",
            [&parameters](int val) { parameters.set_verbosity_level(val); });
        lua["pic_config"] = pic_config;

        // --- EXECUTE CONFIGURATION SCRIPT ---
        sol::state_view sv(static_cast<lua_State*>(lua));
        auto            result = sv.script_file(PIC::Config::cfg_script_name());

        if (!result.valid()) {
            sol::error  err       = result;
            std::string error_msg = fmt::format("FATAL: Lua execution failed: {}\n", err.what());
            fmt::print(stderr, "{}", error_msg);
            throw std::runtime_error(error_msg);
        }

        // --- POST-SCRIPT SYNCHRONIZATION ---
        // Read electron model mode from Lua global (set in main.lua).
        if (lua["cold_electrons_enabled"] != sol::nil) {
            parameters.cold_electrons_enabled = lua["cold_electrons_enabled"].get<bool>();
        }

        parameters.h                   = grid.step();
        parameters.grid_size_x         = grid.size_x();
        parameters.grid_size_y         = grid.size_y();
        parameters.grid_size_z         = grid.size_z();
        parameters.total_particles_num = particles.total_size();

        if (!parameters.is_valid()) {
            std::cerr << parameters;
            throw std::runtime_error("Invalid parameters after Lua configuration.");
        }

        // Finalize logging (applicability checks + lua_params.txt) now that
        // the ENTIRE script has run, including any pic_parameters overrides
        // a case makes after sim_core.run() returns. This guarantees the log
        // always reflects what the simulation actually uses -- see
        // opic_finalize_logging() in sim_core.lua for why this must happen
        // here rather than inside run() itself.
        if (lua["opic_finalize_logging"] != sol::nil) {
            lua["opic_finalize_logging"]();
        }
    } catch (const std::exception& e) {
        std::cerr << "\nStandard Exception: " << e.what() << std::endl;
        return false;
    } catch (...) {
        std::cerr << "\nUndefined error occurred in bind_to_lua." << std::endl;
        return false;
    }

    return true;
}
