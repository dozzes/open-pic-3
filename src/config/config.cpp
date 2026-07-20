#include "config/config.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <stdexcept>
#include <functional>

// static
PIC::Config::Parameters PIC::Config::params = Config::Parameters();

namespace {

static std::ofstream out("opic_trace.log");

// All setters below guard strictly-positive quantities (steps, densities,
// scales): zero would silently produce division by zero downstream
// (inv_h, tau_2h, output normalization), so it is rejected like negatives.
void check_if_positive(const std::string& msg, double val)
{
    if (val <= 0.0)
        throw std::invalid_argument(msg);
}

} // unnamed namespace

// static
std::ofstream& PIC::Config::ofs_log = out;

void PIC::Config::set_dens_cutoff(double dens_cutoff)
{
    check_if_positive("min NP must be > 0.0.", dens_cutoff);
    params.dens_cutoff = dens_cutoff;
}

void PIC::Config::set_tau(double tau)
{
    check_if_positive("tau must be > 0.0.", tau);
    params.tau = tau;
}

void PIC::Config::set_L_scale(double L_scale)
{
    check_if_positive("L_scale must be > 0.0.", L_scale);
    params.L_scale = L_scale;
}

void PIC::Config::set_T_scale(double T_scale)
{
    check_if_positive("T_scale must be > 0.0.", T_scale);
    params.T_scale = T_scale;
}

void PIC::Config::set_U_scale(double U_scale)
{
    check_if_positive("U_scale must be > 0.0.", U_scale);
    params.U_scale = U_scale;
}

void PIC::Config::set_N_scale(double N_scale)
{
    check_if_positive("N_scale must be > 0.0.", N_scale);
    params.N_scale = N_scale;
}

void PIC::Config::set_E_scale(double E_scale)
{
    check_if_positive("E_scale must be > 0.0.", E_scale);
    params.E_scale = E_scale;
}

void PIC::Config::set_B_scale(double B_scale)
{
    check_if_positive("B_scale must be > 0.0.", B_scale);
    params.B_scale = B_scale;
}

// static
void PIC::Config::to_stream(std::ostream& os)
{
    os << params;
}

PIC::Config::Parameters::Parameters()
{
    set_os_name();
}

bool PIC::Config::Parameters::is_valid() const
{
    const bool is_density_threshold_valid = (grid_threshold == GridThreshold::Min_Density) ? (dens_cutoff > 0.0) : true;

    bool valid = (h > 0.0 && tau > 0.0 && is_density_threshold_valid && time_steps > 0 && save_time_steps > 0
                  && grid_size_x > 0 && grid_size_y > 0 && grid_size_z > 0 && L_scale > 0.0 && T_scale > 0.0 && U_scale > 0.0
                  && N_scale > 0.0 && E_scale > 0.0 && B_scale > 0.0);

    // Parameter conflict diagnostics (warnings, not errors).
    if (Spitzer_eta_floor_mult > 1e-3 && cold_electrons_enabled) {
        std::cerr << "\nWARNING: Spitzer_eta_floor_mult=" << Spitzer_eta_floor_mult
                  << " is set, but cold_electrons_enabled=true.\n"
                  << "         FRONTFLOOR has NO EFFECT in COLD mode (no Spitzer eta base).\n"
                  << "         Either: (1) set cold_electrons_enabled=false for TE mode,\n"
                  << "            or (2) remove Spitzer_eta_floor_mult if using COLD.\n\n";
    }

    if (isotropic_curl_enabled && h > 3.0 * L_scale) {
        std::cerr << "\nWARNING: isotropic_curl_enabled=true but h/d_i=" << (h / L_scale)
                  << " >> 1.\n"
                  << "         Benefit of isotropic curl minimal at coarse grids.\n\n";
    }

    return valid;
}

bool PIC::Config::is_on_save_step()
{
    return ((current_time_step() % save_time_steps() == 0) || (current_time_step() == time_steps()));
}

void PIC::Config::set_os_name()
{
#ifdef _WIN64
    params.os_name = "Windows 64-bit";
#elif _WIN32
    params.os_name = "Windows 32-bit";
#elif __unix || __unix__
    params.os_name = "Unix";
#elif __APPLE__ || __MACH__
    params.os_name = "Mac OSX";
#elif __linux__
    params.os_name = "Linux";
#elif __FreeBSD__
    params.os_name = "FreeBSD";
#else
    params.os_name = "Unknown";
#endif
}

namespace PIC {

std::ostream& operator<<(std::ostream& os, const Config::Parameters& params)
{
    os << "\nOpenPIC runtime configuration"
       << "\n----------------------------------------------------------------"
       << "\n  script                 " << params.cfg_script_name
       << "\n  grid                   " << params.grid_size_x << " x " << params.grid_size_y << " x " << params.grid_size_z
       << "\n  h                      " << params.h
       << "\n  tau                    " << params.tau
       << "\n  time steps             " << params.time_steps
       << "\n  save interval          " << params.save_time_steps
       << "\n  total particles        " << params.total_particles_num
       << "\n  density cutoff         " << params.dens_cutoff
       << "\n  process                " << params.process_rank << " / " << params.process_num
       << "\n  OS                     " << params.os_name << std::boolalpha
       << "\n\nDiagnostics"
       << "\n  save particles         " << params.save_particle_diagnostics
       << "\n  save all particles     " << params.save_all_particles
       << "\n  save whole grid        " << params.save_whole_grid
       << "\n  save X planes          " << params.save_grid_x_plains
       << "\n  save Y planes          " << params.save_grid_y_plains
       << "\n  save Z planes          " << params.save_grid_z_plains
       << "\n  raw grid debug         " << params.save_raw_grid_debug
       << "\n\nScales"
       << "\n  L_scale                " << params.L_scale
       << "\n  T_scale                " << params.T_scale
       << "\n  U_scale                " << params.U_scale
       << "\n  N_scale                " << params.N_scale
       << "\n  E_scale                " << params.E_scale
       << "\n  B_scale                " << params.B_scale
       << "\n\nNumerics"
       << "\n  particle pusher        "
       << (params.particle_push_alg == ParticlePushAlg::Direct  ? "Direct"
           : params.particle_push_alg == ParticlePushAlg::Boris ? "Boris"
                                                                : "Undefined")
       << "\n  scatter                "
       << (params.scatter_alg == ScatterAlg::Standard    ? "Standard"
           : params.scatter_alg == ScatterAlg::Zigzag    ? "Zigzag"
           : params.scatter_alg == ScatterAlg::Esirkepov ? "Esirkepov"
           : params.scatter_alg == ScatterAlg::NGP       ? "NGP"
           : params.scatter_alg == ScatterAlg::TSC       ? "TSC"
                                                         : "Undefined")
       << "\n  grid threshold         "
       << (params.grid_threshold == GridThreshold::Min_Density ? "Min_Density"
           : params.grid_threshold == GridThreshold::Local_CFL ? "Local_CFL"
                                                               : "Undefined")
       << "\n  CFL severity           "
       << (params.CFL_severity == CFLSeverity::Ignore   ? "Ignore"
           : params.CFL_severity == CFLSeverity::Absorb ? "Absorb"
           : params.CFL_severity == CFLSeverity::Stop   ? "Stop"
                                                        : "Undefined")
       << "\n  magnetic solver        "
       << (params.magnetic_field_alg == MagneticFieldAlg::FDTD   ? "FDTD"
           : params.magnetic_field_alg == MagneticFieldAlg::PSTD ? "PSTD"
                                                                 : "Undefined")
       << "\n  random seed            " << params.random_seed
       << "\n  isotropic curl         " << params.isotropic_curl_enabled
       << "\n  cold electrons         " << params.cold_electrons_enabled
       << "\n  filtering              " << params.use_filtering
       << "\n  electron heat chi      " << params.electron_thermal_conductivity
       << "\n  resistivity            " << params.resistivity
       << "\n  Spitzer Te ref         " << params.Spitzer_Te_ref
       << "\n  Spitzer Te smooth pass " << params.Spitzer_Te_smooth_passes
       << "\n  Spitzer eta floor mult " << params.Spitzer_eta_floor_mult
       << "\n  hyper resistivity      " << params.hyper_resistivity
       << "\n  kappa friction         " << params.kappa_friction
       << "\n  kappa B ref            " << params.kappa_B_ref
       << "\n\nConstants"
       << "\n  c                      " << Constants::c()
       << "\n  e                      " << Constants::e()
       << "\n  mp                     " << Constants::mp()
       << "\n  e/mp                   " << Constants::e_mp()
       << "\n  pi                     " << Constants::pi()
       << "\n  c/4pi/e                " << Constants::c_4pi_e()
       << "\n----------------------------------------------------------------"
       << std::endl;

    return os;
}

} // namespace PIC
