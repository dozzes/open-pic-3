#pragma once

#include "util/opic_fwd.h"
#include "config/constants.h"

#include <string>
#include <iosfwd>

namespace PIC {

enum class ParticlePushAlg { Direct = 0, Boris };
// Selects the particle shape used consistently for particle-to-grid
// deposition and grid-to-particle field gathering.  `Standard` is the
// historical CIC implementation; NGP and quadratic TSC are controls for
// particle-shape and hybrid particle/fluid cancellation errors.
enum class ScatterAlg { Standard = 0, Zigzag, Esirkepov, NGP, TSC };
enum class GridThreshold { Min_Density = 0, Local_CFL };
enum class CFLSeverity { Ignore, Absorb, Stop };
enum class MagneticFieldAlg { FDTD = 0, PSTD };

class Config
{
  public:
    struct Parameters
    {
        Parameters();

        bool is_valid() const;

        // Derived step coefficients: the single source for solver prefactors.
        // Member functions (not the Config static wrappers) so that solver
        // code taking `const Parameters&` explicitly needs no global state.
        double inv_h() const noexcept { return 1.0 / h; }
        double h_2() const noexcept { return 0.5 * h; }
        double tau_2() const noexcept { return 0.5 * tau; }
        double tau_2h() const noexcept { return tau_2() / h; }
        double ctau_2() const noexcept { return 0.5 * Constants::c() * tau; }
        double ctau_h() const noexcept { return Constants::c() * tau / h; }
        double ctau_2h() const noexcept { return 0.5 * ctau_h(); }
        double c_4pi_e_h() const noexcept { return Constants::c_4pi_e() / h; }

        std::string cfg_script_name;
        double  h = 0.0;               // spatial step
        double  tau = 0.0;             // time step
        double  dens_cutoff = 0.0;     // minimum possible value of grid NP (density cutoff value)
        bool    use_filtering = false;
        index_t time_steps = 0;      // time steps number
        index_t save_time_steps = 0; // save time steps number

        index_t grid_size_x = 0;
        index_t grid_size_y = 0;
        index_t grid_size_z = 0;

        index_t total_particles_num = 0;

        unsigned long long random_seed = 42; // Seed for the mt19937_64 RNG (utility.h).
                                        // Fixed default keeps runs reproducible by default;
                                        // must be set (from Lua) before the first rnd() call
                                        // to take effect -- engines seed lazily per thread.

        bool save_all_particles = false;
        bool save_particle_diagnostics = true;
        bool save_raw_grid_debug = false;
        bool save_whole_grid = false;
        bool save_grid_x_plains = false;
        bool save_grid_y_plains = false;
        bool save_grid_z_plains = false;

        std::string os_name;
        index_t process_num = 1;
        index_t process_rank = 0;

        index_t current_time_step = 0;

        double L_scale = 1.0;
        double T_scale = 1.0;
        double U_scale = 1.0;
        double N_scale = 1.0;
        double E_scale = 1.0;
        double B_scale = 1.0;

        ParticlePushAlg particle_push_alg = ParticlePushAlg::Direct;
        ScatterAlg      scatter_alg = ScatterAlg::Standard;
        GridThreshold   grid_threshold = GridThreshold::Min_Density;
        double          electron_thermal_conductivity = 0.0;
        double          resistivity = 0.0;
        double          Spitzer_Te_ref = 0.0; // Reference temperature (K) at which resistivity is defined.
                                        // eta(Te) = resistivity * (Spitzer_Te_ref / Te)^1.5
                                        // Set to 0 to disable Spitzer scaling (constant resistivity).
        double          hyper_resistivity = 0.0; // eta4 [s*cm^2]: induction filter
                                        // dB/dt = -(eta4*c^2/4pi)*Lap^2(B), applied after the
                                        // B corrector (apply_hyper_resistivity). Damping ~ k^4
                                        // selectively removes grid-scale (2h) modes (m=4 seed,
                                        // checkerboard). NOT a term in Ohm's law E: pushing
                                        // particles with it injects unphysical forces where shot
                                        // noise reseeds grid modes. 0 disables the filter.
        int             Spitzer_Te_smooth_passes = 0; // Number of 3-tap binomial
                                        // ([1,2,1]/4) smoothing passes applied to Te before it
                                        // enters eta(Te) in update_Spitzer_coefficients. A single
                                        // pass exactly cancels a period-2h checkerboard in Te
                                        // while barely touching structure spanning several cells
                                        // (e.g. the compressed front). The real cell.Te is left
                                        // untouched -- only the value fed to the eta formula is
                                        // filtered. 0 disables filtering (eta uses raw cell.Te,
                                        // as before); only takes effect when Spitzer_Te_ref > 0.
        double          Spitzer_eta_floor_mult = 1e-3; // eta_floor = base_eta * this.
                                        // Where Te > Spitzer_Te_ref (compression heating at the
                                        // front/shell) the adaptive formula gives eta BELOW
                                        // base_eta -- measured at Ma=20.3 to be as low as ~0.5x
                                        // base_eta right in the shell where the m=4 mode grows,
                                        // while the cold cavity core independently sits at the
                                        // opposite extreme (eta_ceil, ~2.5x base_eta). Raising
                                        // this to 1.0 pins the front/shell's eta at base_eta
                                        // (like Spitzer_Te_ref=0) while leaving the cold cavity's
                                        // physically-adaptive ceiling untouched -- a targeted
                                        // alternative to freezing eta(Te) everywhere.
                                        // Default 1e-3 reproduces the original permissive floor.
        double          kappa_friction = 0.0; // Dimensionless ion-electron friction
                                        // coefficient, legacy-Fortran-style ("xappa"). When > 0
                                        // (and kappa_B_ref > 0) it REPLACES the Spitzer eta(Te)
                                        // model: E_res = kappa*(B_ref/c)*(Ui-Ue), implemented as
                                        // a per-cell resistivity
                                        //     eta = kappa * kappa_B_ref / (c * e * ne)
                                        // (since Ui-Ue = j/(e*ne) and E_res = eta*j). Unlike the
                                        // Spitzer model this eta has NO Te feedback and scales as
                                        // 1/ne: weak in the dense shell/core, strong in rarefied
                                        // regions -- the legacy code's fielde divides its resistive term
                                        // curl(H)/n by local density. The same cell.eta feeds the
                                        // Joule term of apply_electron_thermodynamics, which then
                                        // reproduces the legacy code's temp2 friction heating
                                        // 0.5*tau*xappa*|Ui-Ue|^2. Additionally every ion feels
                                        // the reaction drag du/dt = kappa*Omega_ci*(Ue-u)
                                        // (apply_kappa_ion_friction, the move_ions channel) --
                                        // the ion-side half of the same interspecies friction,
                                        // which standard hybrid formulations drop.
                                        // 0 disables all three channels (default).
        double          kappa_B_ref = 0.0; // [G] Reference field B0 for kappa_friction --
                                        // the legacy code's field unit h0 (the uniform background field,
                                        // i.e. Bz0 for this problem). Must be > 0 for the kappa
                                        // mode to engage.
        bool            isotropic_curl_enabled = false; // Transversely-smoothed
                                        // (Cole-Karkkainen-style) curl stencils in the field
                                        // solver: every first difference D_a is replaced by
                                        // D_a*(1 + g*(d2_b + d2_c)), g = 1/24, where d2 are the
                                        // second differences in the two transverse axes. Applied
                                        // consistently to curl(E) (Faraday) and curl(B) (Ampere Ue
                                        // and the resistive Ohm term), it cancels the leading
                                        // O(h^2*k^4) ANISOTROPIC term of the discrete curl-curl
                                        // operator -- the suspected driver of the m=4 square mode
                                        // on sharp fronts -- while keeping div(curl)=0 exactly
                                        // (the smoothing factors commute with the differences).
                                        // The isotropic part of the h^2 dispersion error remains.
                                        // Default false reproduces the classic 2-point Yee curl
                                        // bit-for-bit.
        bool            cold_electrons_enabled = true; // Electron model: cold (no resistivity/heat)
                                        // vs. thermodynamic (Spitzer resistivity, kappa conduction).
        CFLSeverity      CFL_severity = CFLSeverity::Stop;
        MagneticFieldAlg magnetic_field_alg = MagneticFieldAlg::FDTD;
        int             verbosity_level = 1; // Profiling/diagnostic output level:
                                        // 0 = silent (no timing output)
                                        // 1 = throttled progress with ETA and throughput
                                        // 2 = progress + intermediate timers (detailed profiling)

        // Setter methods for Lua configuration (called from main.lua before sim_core).
        // Each parameter can be set individually without touching others.
        void set_Spitzer_eta_floor_mult(double val) { Spitzer_eta_floor_mult = val; }
        void set_Spitzer_Te_ref(double val) { Spitzer_Te_ref = val; }
        void set_isotropic_curl_enabled(bool val) { isotropic_curl_enabled = val; }
        void set_cold_electrons_enabled(bool val) { cold_electrons_enabled = val; }
        void set_kappa_friction(double val) { kappa_friction = val; }
        void set_kappa_B_ref(double val) { kappa_B_ref = val; }
        void set_resistivity(double val) { resistivity = val; }
        void set_verbosity_level(int val) { verbosity_level = val; }
    };

    static const std::string& cfg_script_name() { return params.cfg_script_name; }

    static void   set_tau(double tau);
    static double tau() noexcept { return params.tau; }
    static double h() noexcept { return params.h; }
    static double inv_h() noexcept { return params.inv_h(); }
    static double h_2() noexcept { return params.h_2(); }
    static void   set_dens_cutoff(double dens_cutoff); // min density (NP)
    static double dens_cutoff() { return params.dens_cutoff; }
    static void   set_use_filtering(bool use_filtering) { params.use_filtering = use_filtering; }
    static bool   use_filtering() { return params.use_filtering; }

    static double tau_2() noexcept { return params.tau_2(); }
    static double tau_2h() noexcept { return params.tau_2h(); }
    static double ctau_2() noexcept { return params.ctau_2(); }
    static double ctau_h() noexcept { return params.ctau_h(); }
    static double ctau_2h() noexcept { return params.ctau_2h(); }
    static double c_4pi_e_h() noexcept { return params.c_4pi_e_h(); }

    static index_t grid_size_x() noexcept { return params.grid_size_x; }
    static index_t grid_size_y() noexcept { return params.grid_size_y; }
    static index_t grid_size_z() noexcept { return params.grid_size_z; }

    static index_t total_particles_num() { return params.total_particles_num; }

    static unsigned long long random_seed() noexcept { return params.random_seed; }

    static bool save_all_particles() { return params.save_all_particles; }
    static bool save_particle_diagnostics() { return params.save_particle_diagnostics; }
    static bool save_raw_grid_debug() { return params.save_raw_grid_debug; }
    static bool save_whole_grid() { return params.save_whole_grid; }
    static bool save_grid_x_plains() { return params.save_grid_x_plains; }
    static bool save_grid_y_plains() { return params.save_grid_y_plains; }
    static bool save_grid_z_plains() { return params.save_grid_z_plains; }

    static std::string os_name() { return params.os_name; }
    static index_t process_num() { return params.process_num; }
    static index_t process_rank() { return params.process_rank; }

    static index_t time_steps() { return params.time_steps; }           // time steps number
    static index_t save_time_steps() { return params.save_time_steps; } // save time steps number
    static index_t current_time_step() { return params.current_time_step; }
    static void    set_current_time_step(index_t time_step) { params.current_time_step = time_step; }
    static bool    is_on_save_step();
    static void    set_os_name();

    static double L_scale() { return params.L_scale; }
    static double T_scale() { return params.T_scale; }
    static double U_scale() { return params.U_scale; }
    static double N_scale() { return params.N_scale; }
    static double E_scale() { return params.E_scale; }
    static double B_scale() { return params.B_scale; }

    static void set_L_scale(double L_scale);
    static void set_T_scale(double T_scale);
    static void set_U_scale(double U_scale);
    static void set_N_scale(double N_scale);
    static void set_E_scale(double E_scale);
    static void set_B_scale(double B_scale);

    static void            set_particle_push_alg(ParticlePushAlg alg) { params.particle_push_alg = alg; }
    static ParticlePushAlg particle_push_alg() { return params.particle_push_alg; }

    static void       set_scatter_alg(ScatterAlg alg) { params.scatter_alg = alg; }
    static ScatterAlg scatter_alg() { return params.scatter_alg; }

    static void          set_grid_threshold(GridThreshold gt) { params.grid_threshold = gt; }
    static GridThreshold grid_threshold() { return params.grid_threshold; }

    static void        set_CFL_severity(CFLSeverity CFL_severity) { params.CFL_severity = CFL_severity; }
    static CFLSeverity CFL_severity() { return params.CFL_severity; }

    static void             set_magnetic_field_alg(MagneticFieldAlg alg) { params.magnetic_field_alg = alg; }
    static MagneticFieldAlg magnetic_field_alg() { return params.magnetic_field_alg; }

    static bool isotropic_curl_enabled() { return params.isotropic_curl_enabled; }
    static bool cold_electrons_enabled() { return params.cold_electrons_enabled; }

    // No setter: Lua writes Parameters::electron_thermal_conductivity directly
    // (bind_to_lua), and 0.0 is a legal value (disables the heat solver).
    static double electron_thermal_conductivity() noexcept { return params.electron_thermal_conductivity; }

    static double hyper_resistivity() noexcept { return params.hyper_resistivity; }

    static Parameters& parameters() { return params; }

    static void to_stream(std::ostream& os);

    static std::ofstream& logger() { return ofs_log; }

  private:
    static Parameters     params;
    static std::ofstream& ofs_log;
};

std::ostream& operator<<(std::ostream& out, const Config::Parameters& params);

} // namespace PIC
