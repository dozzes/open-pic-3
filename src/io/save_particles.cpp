#include "io/save_particles.h"
#include "grid/grid.h"
#include "particles/particles.h"
#include "core/gather_scatter.h"
#include "particles/particle_groups.h"
#include "particles/marker_particles.h"
#include "core/gather_scatter.h"
#include "config/config.h"
#include "io/io_utilities.h"

#include <omp.h>
#include <cstdio>

using namespace std;

namespace {

// snprintf into a caller-owned buffer instead of fprintf: lets the per-particle
// formatting run in parallel (each thread formats into its own string), with a
// single sequential fwrite per thread buffer at the end doing all the actual I/O.
int format_particle_line(const Grid& grid, const Particles& particles, index_t p, char* buf, size_t buf_size)
{
    Grid::NodeType  point_val;
    const Particle& particle = particles[p];
    PIC::from_grid_to_point(grid, particle.r, point_val);

    const double L = PIC::Config::L_scale();
    const double U = PIC::Config::U_scale();
    const double E = PIC::Config::E_scale();
    const double B = PIC::Config::B_scale();

    return snprintf(buf,
                     buf_size,
                     "%e\t%e\t%e\t"  // (particle.r.x/L), particle.r.y/L), (particle.r.z/L)
                     "%e\t"          // (particle.v.abs()/U)
                     "%e\t%e\t%e\t"  // (particle.v.x/U), particle.v.y/U), (particle.r.z/U)
                     "%e\t"          // (point_val.B.abs()/B)
                     "%e\t%e\t%e\t"  // (point_val.B.x/B), (point_val.B.y/B), (point_val.B.z/B)
                     "%e\t"          // (point_val.E.abs()/E)
                     "%e\t%e\t%e\n", // (point_val.B.x/B), (point_val.B.y/B), (point_val.B.z/B)
                     (particle.r.x / L),
                     (particle.r.y / L),
                     (particle.r.z / L),
                     (particle.v.abs() / U),
                     (particle.v.x / U),
                     (particle.v.y / U),
                     (particle.v.z / U),
                     (point_val.B.abs() / B),
                     (point_val.B.x / B),
                     (point_val.B.y / B),
                     (point_val.B.z / B),
                     (point_val.E.abs() / E),
                     (point_val.E.x / E),
                     (point_val.E.y / E),
                     (point_val.E.z / E));
}

void save_particle(const Grid& grid, const Particles& particles, index_t p, FILE* fout)
{
    if (particles.is_inactive(p))
        return;

    char      line[320];
    const int len = format_particle_line(grid, particles, p, line, sizeof(line));
    fwrite(line, 1, len, fout);
}

void save_particles_group(const Grid& grid, const Particles& particles, const string& group_name)
{
    const string file_name = create_out_file_name(group_name, "parts", PIC::Config::current_time_step());

    const int  group_id = ParticleGroups::get_id_by_name(group_name);
    const bool move_all = (group_name == ParticleGroups::all_particles_name);

    const long     n           = (long)particles.size();
    const int      max_threads = omp_get_max_threads();
    vector<string> thread_bufs(max_threads);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        string&   buf = thread_bufs[tid];
        char      line[320];

#pragma omp for schedule(static)
        for (long p = 0; p < n; ++p) {
            const index_t idx = (index_t)p;
            if (particles.is_inactive(idx))
                continue;

            const Particle& particle = particles[idx];
            if (!(move_all || group_id == particle.group_id))
                continue;

            const int len = format_particle_line(grid, particles, idx, line, sizeof(line));
            buf.append(line, len);
        }
    }

    FILE* fout = fopen(file_name.c_str(), "w");

    if (!fout) {
        cerr << "\nCan't save data: " << file_name << endl;
        return;
    }

    for (const string& buf : thread_bufs)
        fwrite(buf.data(), 1, buf.size(), fout);

    fclose(fout);
}
} //  unnamed namespace

void save_particles(const Grid& grid, const Particles& particles)
{
    print_tm("save_particles");

    const vector<string> group_names = ParticleGroups::group_names();

    for (index_t i = 0; i != group_names.size(); ++i) {
        const string&                        group_name = group_names[i];
        const ParticleGroups::ParticleGroup& group      = ParticleGroups::get_group(group_name);
        if (group.diag & PIC::save_positions)
            save_particles_group(grid, particles, group_name);
    }

    MarkerParticles                                   markers;
    MarkerParticles::MapNameToMarkers::const_iterator it = markers.begin();

    for (; it != markers.end(); ++it) {
        const deque<index_t>& vec_idx = (*it).second;

        if (vec_idx.empty())
            continue;

        const string file_name = create_out_file_name((*it).first, "markers", PIC::Config::current_time_step());
        FILE*        fout      = fopen(file_name.c_str(), "w");

        if (!fout) {
            cerr << "\nCan't save data: " << file_name << endl;
        } else {
            const char* the_header = "X\tY\tZ\t"
                                     "V\tVx\tVy\tVz\t"
                                     "B\tBx\tBy\tBz\t"
                                     "E\tEx\tEy\tEz\n";
            fputs(the_header, fout);
            for (index_t i = 0; i != vec_idx.size(); ++i) {
                save_particle(grid, particles, vec_idx[i], fout);
            }
            fclose(fout);
        }
    }

    if (PIC::Config::save_all_particles()) {
        save_particles_group(grid, particles, ParticleGroups::all_particles_name);
    }

} // unnamed namespace
