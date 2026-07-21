#include "core/simulate.h"
#include "config/config.h"
#include "lua/call_lua_function.h"
#include "particles/particles.h"
#include "core/heat_solver.h"
#include "core/check_particle.h"
#include "particles/particle_groups.h"
#include "core/gather_scatter.h"
#include "io/io_utilities.h"
#include "io/save_grid.h"
#include "io/save_particles.h"
#include "mpi/mpi_runtime.h"

#include "lua/use_lua.h"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "pocketfft_hdronly.h"

#include <fmt/core.h>
#include <omp.h>
#include <stdexcept>

#include "core/electron_model.h"
#include "core/filters.h"

namespace PIC {

void push_particle_std(const Grid& grid, Particle& rp)
{
    DblVector Bp;
    gather_edge(grid, rp.r, &Cell::B, Bp);

    DblVector Ep;
    gather_face(grid, rp.r, &Cell::E, Ep);

    const auto&                          id_table   = ParticleGroups::groups_by_id_table();
    const ParticleGroups::ParticleGroup& part_group = *id_table[rp.group_id];

    const double qm = part_group.charge / part_group.mass * Constants::e_mp();

    const double tau_qm   = Config::tau() * qm;
    const double tau_2qmc = 0.5 * tau_qm / Constants::c();

    const double a = tau_2qmc * Bp.x;
    const double b = tau_2qmc * Bp.y;
    const double c = tau_2qmc * Bp.z;

    const double A = rp.v.x + tau_qm * Ep.x + c * rp.v.y - b * rp.v.z;
    const double B = rp.v.y + tau_qm * Ep.y + a * rp.v.z - c * rp.v.x;
    const double C = rp.v.z + tau_qm * Ep.z + b * rp.v.x - a * rp.v.y;

    const double a2 = a * a;
    const double b2 = b * b;
    const double c2 = c * c;

    const double D = 1.0 / (a2 + b2 + c2 + 1.0);

    const double ab = a * b;
    const double ac = a * c;
    const double bc = b * c;

    /*
    Mathematica solution:
    x -> (A + a^2 A + a b B + B c - b C + a c C)/(1 + a^2 + b^2 + c^2),
    y -> (a A b + B + b^2 B - A c + a C + b c C)/(1 + a^2 + b^2 + c^2),
    z -> (A b - a B + a A c + b B c + C + c^2 C)/(1 + a^2 + b^2 + c^2).
    */

    rp.v.x = D * (A + a2 * A + ab * B + B * c - b * C + ac * C);
    rp.v.y = D * (ab * A + B + b2 * B - A * c + a * C + bc * C);
    rp.v.z = D * (A * b - a * B + ac * A + bc * B + C + c2 * C);
}

void push_particle_Boris(const Grid& grid, Particle& particle)
{
    DblVector B;
    gather_edge(grid, particle.r, &Cell::B, B);

    DblVector E;
    gather_face(grid, particle.r, &Cell::E, E);

    const auto&                          id_table   = ParticleGroups::groups_by_id_table();
    const ParticleGroups::ParticleGroup& part_group = *id_table[particle.group_id];

    const double qm       = part_group.charge / part_group.mass * Constants::e_mp();
    const double tau_2qm  = 0.5 * Config::tau() * qm;
    const double tau_2qmc = tau_2qm / Constants::c();

    const DblVector Vm = particle.v + tau_2qm * E;
    const DblVector V0 = Vm + tau_2qmc * (Vm % B);
    const double    d  = 2.0 * tau_2qmc / (1.0 + tau_2qmc * tau_2qmc * B.get_sqr_len());
    const DblVector Vp = Vm + d * (V0 % B);

    particle.v = Vp + tau_2qm * E;
}

// Ion-side channel of the kappa_friction closure (legacy_code move_ions,
// lines 554-663; see update_Spitzer_coefficients for the E-field channel).
// Interspecies friction acts on BOTH fluids: the +kappa*(Ui-Ue) term in
// Ohm's law is the force on electrons, and momentum conservation requires
// each ion to feel the reaction drag toward the electron fluid:
//     du/dt = nu * (Ue - u),   nu = kappa * Omega_ci = kappa*(q/m)*B_ref/c
// (the legacy code's time unit is t0 = 1/Omega_ci, so its dimensionless
// tau*xappa is exactly tau*nu here; physically nu ~ (m_e/m_i)*nu_ei, the
// ion-electron momentum-exchange rate). Standard hybrid codes -- OpenPIC's
// default included -- drop this force; the legacy code keeps it, applied implicitly per
// particle with denominator s7 = 1 + tau*xappa. We reproduce that as an
// exact implicit sub-step after the E+B push (operator splitting, differs
// from the monolithic solve only at O(tau^2*nu*Omega)): unconditionally
// stable, no new CFL constraint.
void apply_kappa_ion_friction(const Grid& grid, Particle& particle)
{
    const double kappa   = Config::parameters().kappa_friction;
    const double kappa_B = Config::parameters().kappa_B_ref;
    if (kappa <= 0.0 || kappa_B <= 0.0)
        return;

    const auto&                          id_table   = ParticleGroups::groups_by_id_table();
    const ParticleGroups::ParticleGroup& part_group = *id_table[particle.group_id];

    const double qm     = part_group.charge / part_group.mass * Constants::e_mp();
    const double tau_nu = Config::tau() * kappa * qm * kappa_B / Constants::c();

    DblVector UEp;
    gather_face(grid, particle.r, &Cell::UE, UEp);

    particle.v = (1.0 / (1.0 + tau_nu)) * (particle.v + tau_nu * UEp);
}

void push_particle(const Grid& grid, Particle& particle)
{
    if (Config::particle_push_alg() == ParticlePushAlg::Boris)
        push_particle_Boris(grid, particle);
    else
        push_particle_std(grid, particle);

    apply_kappa_ion_friction(grid, particle);
}

void scatter_particle(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    if (Config::scatter_alg() == ScatterAlg::NGP) {
        scatter_particle_ngp(particle, grid, dr);
        return;
    }
    if (Config::scatter_alg() == ScatterAlg::TSC) {
        scatter_particle_tsc(particle, grid, dr);
        return;
    }

    // Zigzag and Esirkepov are disabled: they don't accumulate UP_NP,
    // making UP normalization incorrect at boundaries in the hybrid scheme.
    scatter_particle_std(particle, grid, dr);
}

void scatter_particle_snapshot(const Particle& particle, DensityGrid& grid)
{
    // Diagnostics and t=0 initial moments need instantaneous velocity moments.
    // Charge-conserving scatterers use the particle displacement; with dr=0
    // their current contribution is zero by construction.
    if (Config::scatter_alg() == ScatterAlg::NGP)
        scatter_particle_ngp(particle, grid, DblVector(0, 0, 0));
    else if (Config::scatter_alg() == ScatterAlg::TSC)
        scatter_particle_tsc(particle, grid, DblVector(0, 0, 0));
    else
        scatter_particle_std(particle, grid, DblVector(0, 0, 0));
}

namespace {

#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
struct HighLevelProfileStats
{
    std::vector<std::pair<std::string, double>> seconds_by_name;

    void add(const std::string& name, double seconds)
    {
        for (auto& item : seconds_by_name) {
            if (item.first == name) {
                item.second += seconds;
                return;
            }
        }
        seconds_by_name.emplace_back(name, seconds);
    }

    void print() const
    {
        double total = 0.0;
        for (const auto& item : seconds_by_name)
            total += item.second;

        fmt::print("simulation high-level timings\n");
        for (const auto& item : seconds_by_name) {
            const double pct = total > 0.0 ? (100.0 * item.second / total) : 0.0;
            fmt::print("  {}={:.6f} ({:.1f}%)\n", item.first, item.second, pct);
        }
        fmt::print("  measured_total={:.6f}\n", total);
    }
};

struct ParticlePhaseStats
{
    double can_move_seconds = 0.0;
    double push_seconds     = 0.0;
    double check_seconds    = 0.0;
    double boundary_seconds = 0.0;
    double scatter_seconds  = 0.0;

    long iterations = 0;
    long movable    = 0;
    long pushed     = 0;
    long scattered  = 0;
    long removed    = 0;
    long errors     = 0;
};

void print_particle_phase_stats(const std::string& group_name, const std::vector<ParticlePhaseStats>& stats)
{
    ParticlePhaseStats sum;
    ParticlePhaseStats max_thread;
    double max_thread_total_seconds = 0.0;

    for (const auto& s : stats) {
        const double thread_total_seconds =
            s.can_move_seconds + s.push_seconds + s.check_seconds + s.boundary_seconds + s.scatter_seconds;

        sum.can_move_seconds += s.can_move_seconds;
        sum.push_seconds += s.push_seconds;
        sum.check_seconds += s.check_seconds;
        sum.boundary_seconds += s.boundary_seconds;
        sum.scatter_seconds += s.scatter_seconds;
        sum.iterations += s.iterations;
        sum.movable += s.movable;
        sum.pushed += s.pushed;
        sum.scattered += s.scattered;
        sum.removed += s.removed;
        sum.errors += s.errors;

        max_thread.can_move_seconds = std::max(max_thread.can_move_seconds, s.can_move_seconds);
        max_thread.push_seconds = std::max(max_thread.push_seconds, s.push_seconds);
        max_thread.check_seconds = std::max(max_thread.check_seconds, s.check_seconds);
        max_thread.boundary_seconds = std::max(max_thread.boundary_seconds, s.boundary_seconds);
        max_thread.scatter_seconds = std::max(max_thread.scatter_seconds, s.scatter_seconds);
        max_thread_total_seconds = std::max(max_thread_total_seconds, thread_total_seconds);
        max_thread.iterations = std::max(max_thread.iterations, s.iterations);
        max_thread.scattered = std::max(max_thread.scattered, s.scattered);
    }

    const double total_phase_seconds =
        sum.can_move_seconds + sum.push_seconds + sum.check_seconds + sum.boundary_seconds + sum.scatter_seconds;
    const auto pct = [total_phase_seconds](double value) {
        return total_phase_seconds > 0.0 ? (100.0 * value / total_phase_seconds) : 0.0;
    };

    fmt::print("particle phase timings [{}]\n", group_name);
    fmt::print("  iterations={} movable={} pushed={} scattered={} removed={} errors={}\n",
               sum.iterations,
               sum.movable,
               sum.pushed,
               sum.scattered,
               sum.removed,
               sum.errors);
    fmt::print("  summed thread seconds: can_move={:.6f} ({:.1f}%) push={:.6f} ({:.1f}%) check={:.6f} ({:.1f}%) "
               "boundary={:.6f} ({:.1f}%) scatter={:.6f} ({:.1f}%) total={:.6f}\n",
               sum.can_move_seconds,
               pct(sum.can_move_seconds),
               sum.push_seconds,
               pct(sum.push_seconds),
               sum.check_seconds,
               pct(sum.check_seconds),
               sum.boundary_seconds,
               pct(sum.boundary_seconds),
               sum.scatter_seconds,
               pct(sum.scatter_seconds),
               total_phase_seconds);
    fmt::print("  max thread seconds:    can_move={:.6f} push={:.6f} check={:.6f} boundary={:.6f} scatter={:.6f} total={:.6f}; "
               "max_iterations={} max_scattered={}\n",
               max_thread.can_move_seconds,
               max_thread.push_seconds,
               max_thread.check_seconds,
               max_thread.boundary_seconds,
               max_thread.scatter_seconds,
               max_thread_total_seconds,
               max_thread.iterations,
               max_thread.scattered);
}
#endif

void log_error(const char* msg)
{
    const int         tid           = omp_get_thread_num();
    const std::string log_file_name = fmt::format("opic_thread_{}_move_particles_half_time_err.log", tid);
    std::ofstream     ofs_log(log_file_name.c_str(), std::ios_base::app);
    ofs_log << msg << std::endl;
};

} // namespace

// Appends per-group kinetic energies [erg] to energy_kin.txt on save steps.
// W_g = 1/2 * mass_g*mp * sum(ni * |v|^2) with v at the half-time level
// (n+1/2), consistent with the moments scattered in the same loop.
// Serial runs only: in MPI mode particles are rank-distributed and there is
// no scalar reduction helper; per-group diagnostics are limited there anyway.
static void write_group_kinetic_energy(const std::vector<double>& ni_v2_by_group)
{
    const auto& id_table = ParticleGroups::groups_by_id_table();

    FILE* f = fopen("energy_kin.txt", "a");
    if (!f)
        return;

    // In append mode the initial position is implementation-defined until the
    // first write; seek to the end explicitly so the header check is portable.
    fseek(f, 0, SEEK_END);
    if (ftell(f) == 0) {
        fputs("step\tt", f);
        for (const auto* g : id_table)
            fprintf(f, "\tW_%s", g->name.c_str());
        fputc('\n', f);
    }

    // t matches the normalized-time column of energy.txt (step*tau/T_scale),
    // so the two files join on it directly when plotting the energy balance.
    const double T_scale = Config::parameters().T_scale;
    const double t_norm  = (T_scale > 0.0)
                               ? static_cast<double>(Config::current_time_step()) * Config::tau() / T_scale
                               : 0.0;

    fprintf(f, "%lld\t%e", static_cast<long long>(Config::current_time_step()), t_norm);
    for (size_t gi = 0; gi < id_table.size(); ++gi) {
        const double W = 0.5 * id_table[gi]->mass * Constants::mp() * ni_v2_by_group[gi];
        fprintf(f, "\t%e", W);
    }
    fputc('\n', f);
    fclose(f);
}

void move_particles_half_time(std::vector<DensityGrid>& density_grids,
                              const Grid&          grid,
                              Particles&           particles,
                              const std::string&   group_name)
{
    print_tm("move_particles_half_time: ", group_name);
    const size_t inactive_before_cleanup = particles.inactive_count();
    const size_t compacted               = particles.remove_inactives();
    if (inactive_before_cleanup != 0 || compacted != 0) {
        std::cout << "  particles inactive before cleanup [" << group_name << "]: pending=" << inactive_before_cleanup
                  << " remove_inactives=" << compacted << " slots=" << particles.size() << std::endl;
    }
    const size_t inactive_before_move = particles.inactive_count();

    const size_t           particles_num = particles.size();
    const int              group_id      = ParticleGroups::get_id_by_name(group_name);
    const bool             move_all      = (group_name == ParticleGroups::all_particles_name);
    const PIC::CFLSeverity cfl_severity  = Config::CFL_severity();

    // Pre-filter indices so each thread gets equal real work (no skipped iterations).
    // When move_all=true we use the full range directly to avoid the allocation.
    std::vector<index_t> group_indices;
    if (!move_all) {
        group_indices.reserve(particles_num / 8);
        for (index_t p = 0; p < particles_num; ++p)
            if (particles[p].group_id == group_id)
                group_indices.push_back(p);
    }

    const index_t work_size = move_all ? static_cast<index_t>(particles_num) : static_cast<index_t>(group_indices.size());
    long err_count = 0;
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
    std::vector<ParticlePhaseStats> phase_stats(static_cast<size_t>(omp_get_max_threads()));
#endif

    // Kinetic-energy accumulation is enabled only on save steps, so the hot
    // loop pays nothing on ordinary steps (see write_group_kinetic_energy).
    const size_t energy_groups = ParticleGroups::get_groups_count();
    const bool   do_energy     = move_all && Config::is_on_save_step() && !Mpi::enabled() && energy_groups > 0;
    std::vector<std::vector<double>> ke_by_thread;
    if (do_energy)
        ke_by_thread.assign(static_cast<size_t>(omp_get_max_threads()), std::vector<double>(energy_groups, 0.0));

    // Which domain face a removed particle last crossed (x_lo,x_hi,y_lo,y_hi,z_lo,z_hi).
    // One aggregated line per step in opic_trace.log, not one line per particle,
    // so the log stays bounded regardless of how many particles leave at once.
    std::vector<std::array<size_t, 6>> boundary_hits_by_thread(static_cast<size_t>(omp_get_max_threads()),
                                                                std::array<size_t, 6>{});

    // try-catch is per-thread, not per-iteration: the compiler sees a clean loop body
    // and can optimise it freely (no per-iteration unwind frame).
    // nowait drops the implicit barrier from #pragma omp for; the explicit barrier
    // below ensures all threads synchronise before the parallel region closes.
#pragma omp parallel reduction(+ : err_count)
    {
        const int tid = omp_get_thread_num();
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
        ParticlePhaseStats& local_stats = phase_stats[static_cast<size_t>(tid)];
#endif
        try {
#ifdef OPENPIC_PARTICLE_STATIC_SCHEDULE
#pragma omp for schedule(static) nowait
#else
#pragma omp for schedule(guided) nowait
#endif
            for (index_t i = 0; i < work_size; ++i) {
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                ++local_stats.iterations;
#endif
                const index_t p        = move_all ? i : group_indices[i];
                Particle&     particle = particles[p];
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                double t0 = omp_get_wtime();
#endif
                const bool can_move = is_particle_can_move(particles, p, grid);
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                double t1 = omp_get_wtime();
                local_stats.can_move_seconds += t1 - t0;
#endif
                if (can_move) {
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                    ++local_stats.movable;
                    t0 = omp_get_wtime();
#endif
                    push_particle(grid, particle);
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                    t1 = omp_get_wtime();
                    local_stats.push_seconds += t1 - t0;
                    ++local_stats.pushed;
#endif
                    if (do_energy && particle.group_id >= 0
                        && static_cast<size_t>(particle.group_id) < energy_groups)
                        ke_by_thread[static_cast<size_t>(tid)][static_cast<size_t>(particle.group_id)] +=
                            particle.ni * particle.v.get_sqr_len();

                    const DblVector dr = Config::tau_2() * particle.v;

#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                    t0 = omp_get_wtime();
#endif
                    const bool move_ok = check_particle_move(particle, grid, dr);
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                    t1 = omp_get_wtime();
                    local_stats.check_seconds += t1 - t0;
#endif
                    if (!move_ok) {
                        if (cfl_severity == PIC::CFLSeverity::Absorb)
                            particles.remove_later(p);
                        else {
                            ++err_count;
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                            ++local_stats.errors;
#endif
                        }
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                        ++local_stats.removed;
#endif
                    } else {
                        // Particle is known active (passed is_particle_can_move).
                        // Only check new-position boundary — skip redundant is_inactive.
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                        t0 = omp_get_wtime();
#endif
                        const double  h2  = Config::h();
                        const index_t ni2 = to_idx(particle.r.x + dr.x, h2);
                        const index_t nj2 = to_idx(particle.r.y + dr.y, h2);
                        const index_t nk2 = to_idx(particle.r.z + dr.z, h2);
                        const int     kd2 = 1;
                        const bool in_new_cell =
                            ni2 >= kd2 && ni2 < grid.size_x() - kd2 && nj2 >= kd2 && nj2 < grid.size_y() - kd2 && nk2 >= kd2
                            && nk2 < grid.size_z() - kd2 && grid(ni2, nj2, nk2).state() != PIC::cs_absorptive;
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                        t1 = omp_get_wtime();
                        local_stats.boundary_seconds += t1 - t0;
#endif
                        if (in_new_cell) {
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                            t0 = omp_get_wtime();
#endif
                            scatter_particle(particle, density_grids[tid], dr);
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                            t1 = omp_get_wtime();
                            local_stats.scatter_seconds += t1 - t0;
                            ++local_stats.scattered;
#endif
                            particle.r += dr;
                        } else {
                            particles.remove_later(p);
                            // Classify on the raw (double) position, not on ni2/nj2/nk2: those are
                            // size_t, so a particle that crossed the LOW face (negative coordinate)
                            // wraps around to a huge value and would otherwise be miscounted as HIGH.
                            std::array<size_t, 6>& hits = boundary_hits_by_thread[static_cast<size_t>(tid)];
                            const double px = particle.r.x + dr.x;
                            const double py = particle.r.y + dr.y;
                            const double pz = particle.r.z + dr.z;
                            if (px < kd2 * h2) ++hits[0];
                            if (px >= (grid.size_x() - kd2) * h2) ++hits[1];
                            if (py < kd2 * h2) ++hits[2];
                            if (py >= (grid.size_y() - kd2) * h2) ++hits[3];
                            if (pz < kd2 * h2) ++hits[4];
                            if (pz >= (grid.size_z() - kd2) * h2) ++hits[5];
#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
                            ++local_stats.removed;
#endif
                        }
                    }
                }
            }
        } catch (std::exception& e) {
            log_error(e.what());
            ++err_count;
        } catch (...) {
            log_error("Unknown exception in move_particles_half_time");
            ++err_count;
        }
#pragma omp barrier
    }

    const size_t inactive_after_move = particles.inactive_count();
    if (inactive_after_move != inactive_before_move) {
        std::cout << "  particles marked inactive [" << group_name << "]: +" << (inactive_after_move - inactive_before_move)
                  << " pending=" << inactive_after_move << " slots=" << particles.size() << std::endl;
    }

    std::array<size_t, 6> boundary_hits{};
    for (const auto& thread_hits : boundary_hits_by_thread)
        for (int f = 0; f < 6; ++f) boundary_hits[static_cast<size_t>(f)] += thread_hits[static_cast<size_t>(f)];
    if (boundary_hits[0] + boundary_hits[1] + boundary_hits[2] + boundary_hits[3] + boundary_hits[4] + boundary_hits[5] != 0) {
        Config::logger() << "time step = " << Config::current_time_step() << ": particles removed [" << group_name
                         << "]: x_lo=" << boundary_hits[0] << " x_hi=" << boundary_hits[1] << " y_lo=" << boundary_hits[2]
                         << " y_hi=" << boundary_hits[3] << " z_lo=" << boundary_hits[4] << " z_hi=" << boundary_hits[5]
                         << '\n';
    }

#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
    print_particle_phase_stats(group_name, phase_stats);
#endif

    if (do_energy) {
        std::vector<double> ni_v2(energy_groups, 0.0);
        for (const auto& thread_sums : ke_by_thread)
            for (size_t gi = 0; gi < energy_groups; ++gi)
                ni_v2[gi] += thread_sums[gi];
        write_group_kinetic_energy(ni_v2);
    }

    if (err_count != 0)
        throw std::domain_error("Error occurred during the simulation!\nSee opic_thread_*_err.log files for details.");
}

void move_particles_full_time([[maybe_unused]] const Grid& grid, Particles& particles, const std::string& group_name) // step m+1
{
    print_tm("move_particles_full_time: ", group_name);

    const double tau_2         = Config::tau_2();
    const size_t particles_num = particles.size();
    const int    group_id      = ParticleGroups::get_id_by_name(group_name);
    const bool   move_all      = (group_name == ParticleGroups::all_particles_name);

    // Particles that failed the boundary check at the end of move_particles_half_time
    // are already marked inactive. A full grid-cell state lookup here is redundant.
#pragma omp parallel for schedule(static)
    for (index_t p = 0; p < particles_num; ++p) {
        if ((move_all || group_id == particles[p].group_id) && !particles.is_inactive(p)) {
            particles[p].r += tau_2 * particles[p].v;
        }
    }
}

template <class DensityGridType> void set_grid_UP(const Grid& grid, DensityGridType& dg_group)
{
    print_tm("set_grid_UP");

    using DensityNode = typename DensityGridType::NodeType;
    const index_t m   = 1;
    const index_t to_x = grid.size_x() - m;
    const index_t to_y = grid.size_y() - m;
    const index_t to_z = grid.size_z() - m;

    // Each iteration writes only dg_group(i,j,k) and reads (i,j,k) and its
    // -1 neighbors: no cell is written by more than one iteration, so the
    // loop nest is race-free under OpenMP.
#pragma omp parallel for collapse(3)
    for (index_t i = m; i < to_x; ++i)
        for (index_t j = m; j < to_y; ++j)
            for (index_t k = m; k < to_z; ++k) {
                auto& dn = dg_group(i, j, k);

                double NPx = 0.5 * (dn.NP + dg_group(i - 1, j, k).NP);
                double NPy = 0.5 * (dn.NP + dg_group(i, j - 1, k).NP);
                double NPz = 0.5 * (dn.NP + dg_group(i, j, k - 1).NP);

                if constexpr (DensityNode::has_face_density) {
                    // UP components are deposited on faces; use face-deposited charge
                    // as the denominator. The cell-centered NP average is only a
                    // fallback for charge-conserving current scatterers.
                    if (dn.UP_NP.x > 1.0e-6)
                        NPx = dn.UP_NP.x;
                    if (dn.UP_NP.y > 1.0e-6)
                        NPy = dn.UP_NP.y;
                    if (dn.UP_NP.z > 1.0e-6)
                        NPz = dn.UP_NP.z;
                }

                dn.UP.x = (NPx > 1.0e-6) ? (dn.UP.x / NPx) : 0.0;
                dn.UP.y = (NPy > 1.0e-6) ? (dn.UP.y / NPy) : 0.0;
                dn.UP.z = (NPz > 1.0e-6) ? (dn.UP.z / NPz) : 0.0;
            }
}

template <class DensityGridType> void normalize_NP(const Grid& grid, DensityGridType& dg_group)
{
    print_tm("normalize_NP");

    const double  cell_volume = grid.cell_volume();
    const index_t nx          = grid.size_x();
    const index_t ny          = grid.size_y();
    const index_t nz          = grid.size_z();

#pragma omp parallel for collapse(3)
    for (index_t i = 0; i < nx; ++i)
        for (index_t j = 0; j < ny; ++j)
            for (index_t k = 0; k < nz; ++k) {
                dg_group(i, j, k).NP /= cell_volume;
            }
}

void local_Alfven_CFL(Grid& grid, index_t i, index_t j, index_t k)
{
    // No print_tm() here: this runs per-cell inside set_threshold's parallel
    // sweep, and print_tm writes to shared std::cout -- unconditional per-cell
    // calls would both race across threads and spam the console every step.

    Cell& cell = grid(i, j, k);

    const double h = grid.step();
    DblVector    vec_to_point((i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h);
    DblVector    vec_B;
    gather_edge(grid, vec_to_point, &Cell::B, vec_B);
    const double B = vec_B.abs();

    const double mp    = ParticleGroups::max_mp() * PIC::Constants::mp();
    const double pi    = Constants::pi();
    const double tau_2 = Config::tau_2();

    const double min_dens    = 1.5 * B * B * tau_2 * tau_2 / (4 * pi * h * h * mp);
    const double cell_volume = grid.cell_volume();

    // At this stage grid stores NP*cell_volume
    if (cell.NP < min_dens * cell_volume) {
        // '\n' (not endl): this fires per-cell inside set_threshold's grid sweep;
        // flushing every line dominates the sweep on large grids.
        // critical: set_threshold's sweep is now parallelized, and logger() is one
        // shared ofstream -- concurrent operator<< from multiple threads would race.
#pragma omp critical(opic_logger)
        PIC::Config::logger() << "\nlocal_Alfven_CFL:"
                              << "(i,j,k) = (" << i << "," << j << "," << k << ")"
                              << " cell.NP = " << cell.NP << " min_dens = " << min_dens << '\n';

        cell.NP = min_dens * cell_volume;
    }
}

void backgr_fracture(Grid& grid, index_t i, index_t j, index_t k)
{
    Cell&        cell        = grid(i, j, k);
    const double dens_cutoff = Config::dens_cutoff();

    // At this stage grid.NP stores physical density.
    if (cell.NP < dens_cutoff) {
        // '\n' (not endl): fires per-cell, see local_Alfven_CFL above.
        // critical: see local_Alfven_CFL above -- logger() is one shared ofstream.
#pragma omp critical(opic_logger)
        Config::logger() << "time step = " << Config::current_time_step() << ": backgr_fracture: "
                         << "(i,j,k) = (" << i << "," << j << "," << k << "):"
                         << " cell.NP = " << cell.NP << ": dens_cutoff = " << dens_cutoff << '\n';

        cell.NP = dens_cutoff;
    }
}

void grid_threshold(Grid& grid, index_t i, index_t j, index_t k)
{
    if (Config::grid_threshold() == GridThreshold::Local_CFL) {
        local_Alfven_CFL(grid, i, j, k);
    } else {
        backgr_fracture(grid, i, j, k);
    }
}

void set_threshold(Grid& grid)
{
    print_tm("set_threshold");

    const index_t nx = grid.size_x() - 1;
    const index_t ny = grid.size_y() - 1;
    const index_t nz = grid.size_z() - 1;

    // At this stage grid stores physical density. Each cell only reads its own
    // neighbors (gather_edge) and writes its own NP, so this is race-free across
    // (i,j,k) the same way normalize_NP/set_grid_UP are.
#pragma omp parallel for collapse(3) schedule(static)
    for (index_t i = 1; i < nx; ++i)
        for (index_t j = 1; j < ny; ++j)
            for (index_t k = 1; k < nz; ++k) {
                grid_threshold(grid, i, j, k);
            }
}

// --- Helper functions for the simulation loop ---

/**
 * Helper to collect densities from threads into a DensityGrid (for groups)
 */
void collect_to_density_grid(DensityGrid& target, const std::vector<DensityGrid>& densityGrids)
{
    print_tm("collect_to_density_grid");
    target.reset_current();
    target.add_current(densityGrids);
}

/**
 * Helper to collect densities from threads directly into the main Grid
 */
void collect_to_main_grid(Grid& target, const std::vector<DensityGrid>& densityGrids)
{
    print_tm("collect_to_main_grid");
    target.reset_current();
    target.add_current(densityGrids);
}

/**
 * Normalizes ion fluid properties for the main grid
 */
void finalize_main_grid(Grid& grid)
{
    print_tm("finalize_main_grid");

    // 1. Calculate ion velocities (UP) using raw weights.
    // At this stage, grid.NP contains the sum of particle weights (~2.98 for background).
    // Momentum divided by weight gives the correct physical velocity.
    set_grid_UP(grid, grid);
    set_boundary_conditions("on_set_boundary_UP");

    // 2. Normalize NP: convert accumulated weights into physical density.
    // This divides grid.NP by cell_volume, turning 2.98 into 1.0 (or 0.44 into 0.15).
    // After this call, grid.NP is ready for Ohm's law and saving to disk.
    normalize_NP(grid, grid);

    // 3. Set boundary conditions for normalized density.
    set_boundary_conditions("on_set_boundary_NP");

    // 4. Apply density threshold (floor) in physical-density units.
    // Boundary conditions are already applied, so Lua can always work with densities.
    set_threshold(grid);
}

/**
 * Normalizes ion fluid properties for a specific particle group grid
 */
void finalize_group_grid(Grid& grid, DensityGrid& dg_group, const std::string& group_name)
{
    print_tm("finalize_group_grid: " + group_name);

    // 1. Calculate ion velocities (UP) for the group using raw weights.
    set_grid_UP(grid, dg_group);
    set_boundary_conditions("on_set_boundary_group_UP", dg_group, group_name);

    // Note: Group grids do not need thresholding as they are for diagnostics only.

    // 2. Normalize NP: convert accumulated group weights into physical density.
    normalize_NP(grid, dg_group);

    // 3. Set boundary conditions for normalized group density via Lua.
    set_boundary_conditions("on_set_boundary_group_NP", dg_group, group_name);
}

// Scatters instantaneous per-group moments into dg_group and writes each
// group's grid diagnostics. The identical sequence runs at t=0 and on every
// save step. In MPI mode the group moments are summed across ranks before
// finalization (a no-op in serial runs) and only the root rank writes files.
static void save_group_grids(Grid&                     grid,
                             DensityGrid&              dg_group,
                             std::vector<DensityGrid>& per_thread_density_grids,
                             Particles&                particles)
{
    ParticleGroups part_groups;
    for (const std::string& group_name : part_groups.group_names()) {
        if (group_name == ParticleGroups::all_particles_name)
            continue;

        // Clear thread-local buffers
        for (auto& dg : per_thread_density_grids)
            dg.reset_current();

        const int group_id = ParticleGroups::get_id_by_name(group_name);

#pragma omp parallel for
        for (index_t p = 0; p < (index_t)particles.size(); ++p) {
            Particle& particle = particles[p];
            if (particle.group_id == group_id) {
                const int tid = omp_get_thread_num();
                scatter_particle_snapshot(particle, per_thread_density_grids[tid]);
            }
        }

        // Collect group-specific data and finalize its fluid moments
        collect_to_density_grid(dg_group, per_thread_density_grids);
        Mpi::allreduce_density_grid(dg_group);
        finalize_group_grid(grid, dg_group, group_name);

        // Save individual group diagnostics to disk
        if (Mpi::is_root())
            save_grid(grid, dg_group, group_name);
    }
}

/*
Simulation cycle

// 1. Magnetic field predictor (B: n -> n+1/2)
calc_magnetic_field_half_time(grid);
set_boundary_conditions("on_set_boundary_MF");

// 2. Initial resistivity update (eta: n)
update_Spitzer_coefficients(grid);

// 3. Electron Temperature update (Te: n -> n+1)
PIC::HeatSolver::solve(grid);
set_boundary_conditions("on_set_boundary_Te");

// 4. Final resistivity update (eta: n+1)
update_Spitzer_coefficients(grid);

// 5. Electron velocity (Ue: n+1/2)
calc_electrons_velocity(grid);
set_boundary_conditions("on_set_boundary_UE");

// 6. Electric field (Ohm's Law: E n+1/2)
// Now uses updated Te and eta!
calc_electric_field(grid);
set_boundary_conditions("on_set_boundary_EF");

// 7. Particle push (V: n-1/2 -> n+1/2, R: n -> n+1/2)
for (auto& dg : sDensity_grids) dg.reset_current();
move_particles_half_time(sDensity_grids, grid, particles, ParticleGroups::all_particles_name);

// 8. Accumulate data to the main grid: set grid NP and UP arrays.
collect_to_main_grid(grid, sDensity_grids);

// 9: Diagnostic output for specific groups

// 10: Magnetic field corrector (B: n+1/2 -> n+1)
   calc_magnetic_field_half_time(grid);

// 11: Final position update (R: n+1/2 -> n+1): move particles to the final position for the current time step
move_particles_full_time(grid, particles, ParticleGroups::all_particles_name);
*/
namespace {

bool cell_is_finite(const Cell& c)
{
    return std::isfinite(c.NP) && std::isfinite(c.Te) && std::isfinite(c.eta) && is_finite(c.E) && is_finite(c.B)
        && is_finite(c.UE) && is_finite(c.UP);
}

// Zero-tolerance NaN policy: a single non-finite cell invalidates the run, so
// scan every cell instead of probing one hardcoded location (the old
// grid(10,10,10) check missed NaNs anywhere else and only printed a warning).
// Called on save steps only: the O(N^3) sweep is negligible next to the I/O
// done on the same step, and it guarantees corrupted output is never written
// silently -- the run stops with the offending cells named instead.
void check_grid_finite(const Grid& grid)
{
    const index_t nx = grid.size_x();
    const index_t ny = grid.size_y();
    const index_t nz = grid.size_z();

    long long bad_count = 0;
#pragma omp parallel for collapse(3) reduction(+ : bad_count)
    for (index_t i = 0; i < nx; ++i)
        for (index_t j = 0; j < ny; ++j)
            for (index_t k = 0; k < nz; ++k) {
                if (!cell_is_finite(grid(i, j, k)))
                    ++bad_count;
            }

    if (bad_count == 0)
        return;

    // Serial rescan to name the first offenders; cheap since we only get here
    // on a fatal error.
    constexpr long long max_report = 8;
    long long           reported   = 0;
    for (index_t i = 0; i < nx && reported < max_report; ++i)
        for (index_t j = 0; j < ny && reported < max_report; ++j)
            for (index_t k = 0; k < nz && reported < max_report; ++k) {
                const Cell& c = grid(i, j, k);
                if (cell_is_finite(c))
                    continue;
                ++reported;
                log(fmt::format("Non-finite cell ({}, {}, {}): NP={} Te={} eta={} E=({}, {}, {}) "
                                "B=({}, {}, {}) UE=({}, {}, {}) UP=({}, {}, {})",
                                i, j, k, c.NP, c.Te, c.eta, c.E.x, c.E.y, c.E.z, c.B.x, c.B.y, c.B.z,
                                c.UE.x, c.UE.y, c.UE.z, c.UP.x, c.UP.y, c.UP.z),
                    true);
            }

    throw std::domain_error(fmt::format("{} non-finite grid cell(s) detected at step {} (first {} logged).",
                                        bad_count, Config::current_time_step(),
                                        std::min(bad_count, max_report)));
}

} // namespace

void simulate(Grid& grid, Particles& particles)
{
    if (Mpi::is_root())
        std::cout << (Mpi::enabled() ? "\nOpenPIC: MPI replicated-grid simulation started ...\n"
                                     : "\nOpenPIC: Simulation started ...\n");

    if (Mpi::enabled() && Mpi::is_root() && Config::magnetic_field_alg() == MagneticFieldAlg::PSTD) {
        std::cout << "WARNING: PSTD in MPI mode duplicates the full FFT work on every rank.\n";
    }

    // Single explicit parameter set for the whole run: solver functions read
    // it as `const&` instead of reaching for the Config singleton.
    const Config::Parameters& params = Config::parameters();

    // Scratch buffers shared by the solvers; owned here so nothing in the
    // solver layer keeps hidden static state between runs.
    SolverWorkspace workspace;

    const int     max_threads = omp_get_max_threads();
    const index_t steps_num   = params.time_steps;

    omp_set_dynamic(0);
    omp_set_num_threads(max_threads);

#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
    HighLevelProfileStats high_level_profile;
    const auto profile_block = [&high_level_profile](const std::string& name, auto&& fn) {
        const double t0 = omp_get_wtime();
        fn();
        high_level_profile.add(name, omp_get_wtime() - t0);
    };
#else
    const auto profile_block = [](const std::string&, auto&& fn) { fn(); };
#endif

    const auto calc_B_half = [&]() {
        if (params.magnetic_field_alg == MagneticFieldAlg::PSTD)
            calc_magnetic_field_pstd_half_time(grid, params, workspace);
        else
            calc_magnetic_field_half_time(grid, params);
    };

    // 0. Scatter NP and UP at t=0
    print_tm("scatter NP and UP at t=0");

    std::vector<DensityGrid> perThreadDensityGrid(max_threads,
                                                  DensityGrid(grid.size_x(), grid.size_y(), grid.size_z(), grid.step()));
#pragma omp parallel for
    for (index_t p = 0; p < particles.size(); ++p) {
        const int thread_num = omp_get_thread_num();
        // Scatter initial particle positions without displacement
        Particle& particle = particles[p];
        scatter_particle_snapshot(particle, perThreadDensityGrid[thread_num]);
    }

    // Save Step 0 global state
    print_tm("Step 0: Save initial state");
    collect_to_main_grid(grid, perThreadDensityGrid);
    Mpi::allreduce_grid_density(grid); // no-op without MPI
    finalize_main_grid(grid);
    // Particles are rank-distributed in MPI mode; per-particle diagnostics
    // would need a gather step that does not exist, so they are serial-only.
    if (!Mpi::enabled() && Config::save_particle_diagnostics())
        save_particles(grid, particles);
    if (Mpi::is_root())
        save_grid(grid, grid, ParticleGroups::all_particles_name);

    // Grid for the particle group (NP an UP)
    DensityGrid dg_group(grid.size_x(), grid.size_y(), grid.size_z(), grid.step());

    // Save Step 0 group state
    save_group_grids(grid, dg_group, perThreadDensityGrid, particles);

    // Main time-stepping loop
    for (index_t step = 1; step <= steps_num; ++step) {
        Config::set_current_time_step(step);

        // --- STEP 1: Magnetic field predictor (B: n -> n+1/2) ---
        // Faraday's Law update for the first half-step.
        profile_block("lua:on_iteration_begin", [&]() { call_lua_function("on_iteration_begin"); });
        profile_block("magnetic_field_half", [&]() { calc_B_half(); });
        profile_block("boundary:MF", [&]() { set_boundary_conditions("on_set_boundary_MF"); });

        // --- STEP 2: Electron moments and velocity (Ue: n, Pe: n) ---
        // Calculate Ue: n based on CURRENT ion moments and B-field.
        // This ensures that adiabatic terms in thermodynamics use fresh velocity data.
        profile_block("calc_electrons_velocity", [&]() { calc_electrons_velocity(grid, params); });
        profile_block("boundary:UE", [&]() { set_boundary_conditions("on_set_boundary_UE"); });

        // --- STEP 3: Resistivity and electron heat sources (eta: n, Te source terms) ---
        // Ue is fresh here, so adiabatic electron thermodynamics uses current moments.
        profile_block("update_Spitzer_coefficients:first", [&]() { update_Spitzer_coefficients(grid, params, workspace); });
        profile_block("boundary:eta:first", [&]() { set_boundary_conditions("on_set_boundary_eta"); });
        profile_block("apply_electron_thermodynamics", [&]() { apply_electron_thermodynamics(grid, params); });

        // --- STEP 4: Electron temperature diffusion (Te: n -> n+1) ---
        profile_block("heat_solver", [&]() { PIC::HeatSolver::solve(grid, step); });
        profile_block("boundary:Te", [&]() { set_boundary_conditions("on_set_boundary_Te"); });

        // --- STEP 5: Final resistivity update for Ohm's law (eta: n+1) ---
        profile_block("update_Spitzer_coefficients:final", [&]() { update_Spitzer_coefficients(grid, params, workspace); });
        profile_block("boundary:eta:final", [&]() { set_boundary_conditions("on_set_boundary_eta"); });

        // --- STEP 6: Electric field calculation (E: n+1/2)
        // Solve generalized Ohm's law using current magnetic field and electron dynamics.
        profile_block("calc_electric_field", [&]() { calc_electric_field(grid, params); });
        profile_block("boundary:EF", [&]() { set_boundary_conditions("on_set_boundary_EF"); });

        // --- STEP 7: Particle push and scatter (V: n-1/2 -> n+1/2, R: n -> n+1/2) ---
        // Move particles using E: n+1/2 and B: n+1/2.
        // This fills thread-local buffers (density_grids) with data from all particles.
        profile_block("density_grid_reset:main", [&]() {
            for (auto& dg : perThreadDensityGrid)
                dg.reset_current();
        });
        profile_block("move_particles_half_time", [&]() {
            move_particles_half_time(perThreadDensityGrid, grid, particles, ParticleGroups::all_particles_name);
        });

        // --- STEP 8: Accumulate data to the main grid ---
        // Move data from thread-local buffers to the global grid NP and UP arrays.
        // In MPI mode all ranks then sum their density moments (no-op without MPI).
        profile_block("collect_to_main_grid", [&]() { collect_to_main_grid(grid, perThreadDensityGrid); });
        profile_block("mpi:allreduce_density", [&]() { Mpi::allreduce_grid_density(grid); });

        // Post-process grid after all accumulations and filtering
        if (step % 10 == 0 && Config::use_filtering()) {
            profile_block("filter_sources", [&]() { filter_sources(grid); });
        }
        profile_block("finalize_main_grid", [&]() { finalize_main_grid(grid); });

        // Root-only: the grid is replicated across ranks, so one rank's scan
        // covers everything; a root throw reaches Mpi::abort in main().
        if (Config::is_on_save_step() && Mpi::is_root()) {
            profile_block("check_grid_finite", [&]() { check_grid_finite(grid); });
        }

        // --- STEP 9: Diagnostic output for specific groups ---
        if (Config::is_on_save_step()) {
            profile_block("diagnostics:groups",
                          [&]() { save_group_grids(grid, dg_group, perThreadDensityGrid, particles); });
        }

        if (Config::is_on_save_step()) {
            profile_block("save:main", [&]() {
                // Serial-only: see the Step 0 comment on rank-distributed particles.
                if (!Mpi::enabled() && Config::save_particle_diagnostics())
                    save_particles(grid, particles);
                if (Mpi::is_root())
                    save_grid(grid, grid, ParticleGroups::all_particles_name);
            });
        }

        // --- STEP 10: Magnetic field corrector (B: n+1/2 -> n+1) ---
        // Update magnetic field to the full time step
        profile_block("magnetic_field_half", [&]() { calc_B_half(); });
        profile_block("boundary:MF", [&]() { set_boundary_conditions("on_set_boundary_MF"); });
        profile_block("hyper_resistivity", [&]() { apply_hyper_resistivity(grid, params, workspace); });

        // --- STEP 11: Final position update (R: n+1/2 -> n+1) ---
        // Move particles to the final position for the current time step
        profile_block("move_particles_full_time", [&]() {
            move_particles_full_time(grid, particles, ParticleGroups::all_particles_name);
        });
        profile_block("lua:on_particles_moved_full_time", [&]() { call_lua_function("on_particles_moved_full_time"); });

        if (PIC::Config::save_raw_grid_debug() && PIC::Config::current_time_step() < 10 && Mpi::is_root()) {
            profile_block("save_raw_grid", [&]() { save_raw_grid(grid); });
        }

        // Root-only: on_iteration_end typically appends to shared log/diag files,
        // which must not be written by every rank.
        if (Mpi::is_root())
            profile_block("lua:on_iteration_end", [&]() { call_lua_function("on_iteration_end"); });

        if (Mpi::is_root() && should_print_step_marker())
            print_simulation_progress(step, steps_num);
    }

    if (Mpi::is_root())
        finish_simulation_progress();

#ifdef OPENPIC_PROFILE_PARTICLE_PHASES
    if (Mpi::is_root())
        high_level_profile.print();
#endif
}

} // namespace PIC
