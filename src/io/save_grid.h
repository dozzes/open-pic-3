#pragma once

#include "grid/grid_filters.h"
#include "core/gather_scatter.h"
#include "particles/particle_groups.h"
#include "config/constants.h"
#include "io/io_utilities.h"

#include <omp.h>
#include <cstdio>
#include <vector>
#include <string>

template <class GridT, class DensGridT, class PointFilterFoT>
void save_subgrid(const GridT& grid, const DensGridT& dens_grid, const PointFilterFoT& point_filter)
{
    using namespace PIC;

    std::string fileName = create_out_file_name(point_filter.name(), "grd", PIC::Config::current_time_step());

    FILE* fout = fopen(fileName.c_str(), "w");
    if (!fout) {
        std::cerr << "can't save data: " << fileName << std::endl;
        return;
    }

    const char* the_header = "X\tY\tZ\t"
                             "NP\t"
                             "B\tBx\tBy\tBz\t"
                             "E\tEx\tEy\tEz\tdiv(E)\t"
                             "UP\tUPx\tUPy\tUPz\tMa\t"
                             "UE\tUEx\tUEy\tUEz\tTe\tc_Wpi\n";

    fputs(the_header, fout);

    const double h  = grid.step();
    const double L  = Config::L_scale();
    const double U  = Config::U_scale();
    const double N  = Config::N_scale();
    const double E  = Config::E_scale();
    const double B  = Config::B_scale();
    const double Te = Config::parameters().Spitzer_Te_ref;

    const index_t m  = 1;
    const long    nx = (long)(grid.size_x() - m);
    const long    ny = (long)(grid.size_y() - m);
    const long    nz = (long)(grid.size_z() - m);

    // Per-node gather (from_grid_to_point + 6x gather_vector for div(E)) dominates the
    // cost, not the I/O itself: format in parallel into per-thread buffers, then do the
    // actual (sequential, order-preserving) file write once, outside the parallel region.
    const int           max_threads = omp_get_max_threads();
    std::vector<std::string> thread_bufs(max_threads);

#pragma omp parallel
    {
        const int    tid = omp_get_thread_num();
        std::string& buf = thread_bufs[tid];
        char         line[512];

#pragma omp for collapse(3) schedule(static)
        for (long i = (long)m; i < nx; ++i)
            for (long j = (long)m; j < ny; ++j)
                for (long k = (long)m; k < nz; ++k) {
                    DblVector node(i * 1.0, 1.0 * j * 1.0, k * 1.0);
                    DblVector point(i * h, j * h, k * h);

                    if (point_filter(node, point)) {
                        typename GridT::NodeType point_val;
                        from_grid_to_point(grid, point, point_val);

                        typename DensGridT::NodeType dgNode;
                        dgNode.NP = gather_scalar<PIC::CellCentering>(dens_grid, point, &DensGridT::NodeType::NP);

                        gather_face(dens_grid, point, &DensGridT::NodeType::UP, dgNode.UP);

                        double vA = point_val.B.abs() / sqrt(4 * Constants::pi() * point_val.NP * Constants::mp());

                        const double d    = 0.5 * h;
                        const double Ex_r = gather_vector<PIC::FaceXCentering>(
                            grid, point.x + d, point.y, point.z, &GridT::NodeType::E, &DblVector::x);
                        const double Ex_l = gather_vector<PIC::FaceXCentering>(
                            grid, point.x - d, point.y, point.z, &GridT::NodeType::E, &DblVector::x);

                        const double Ey_r = gather_vector<PIC::FaceYCentering>(
                            grid, point.x, point.y + d, point.z, &GridT::NodeType::E, &DblVector::y);
                        const double Ey_l = gather_vector<PIC::FaceYCentering>(
                            grid, point.x, point.y - d, point.z, &GridT::NodeType::E, &DblVector::y);

                        const double Ez_r = gather_vector<PIC::FaceZCentering>(
                            grid, point.x, point.y, point.z + d, &GridT::NodeType::E, &DblVector::z);
                        const double Ez_l = gather_vector<PIC::FaceZCentering>(
                            grid, point.x, point.y, point.z - d, &GridT::NodeType::E, &DblVector::z);

                        const double divE  = (L / E) * ((Ex_r - Ex_l) + (Ey_r - Ey_l) + (Ez_r - Ez_l)) / h;
                        const int    Nz    = 1;
                        const double mi    = 1 * Constants::mp(); //  --background ions mass
                        const double Wpi   = Nz * Constants::e() * 2 * sqrt(Constants::pi() * dgNode.NP / mi);
                        const double c_Wpi = (Wpi > 1e-20) ? (Constants::c() / Wpi) : 0.0;

                        const int len = snprintf(line,
                                sizeof(line),
                                "%e\t%e\t%e\t" /* point.x/L, point.y/L, point.z/L */
                                "%e\t"         /* dgNode.NP/N */
                                "%e\t"         /* point_val.B.abs()/B */
                                "%e\t%e\t%e\t" /* point_val.B.x/B, point_val.B.y/B, point_val.B.z/B */
                                "%e\t"         /* point_val.E.abs()/E */
                                "%e\t%e\t%e\t" /* point_val.E.x/E, point_val.E.y/E, point_val.E.z/E */
                                "%e\t"         /* div(E)/(E/L) */
                                "%e\t"         /* dgNode.UP.abs()/U */
                                "%e\t%e\t%e\t" /* dgNode.UP.x/U, dgNode.UP.y/U, dgNode.UP.z/U */
                                "%e\t"         /* dgNode.UP.abs()/vA */
                                "%e\t"         /* point_val.UE.abs()/U */
                                "%e\t%e\t%e\t" /* dpoint_val.UE.x/U, point_val.UE.y/U, point_val.UE.z/U */
                                "%e\t"         /* point_val.Te */
                                "%e\n",        /* c_Wpi*/
                                (point.x / L),
                                (point.y / L),
                                (point.z / L),
                                (dgNode.NP / N),
                                (point_val.B.abs() / B),
                                (point_val.B.x / B),
                                (point_val.B.y / B),
                                (point_val.B.z / B),
                                (point_val.E.abs() / E),
                                (point_val.E.x / E),
                                (point_val.E.y / E),
                                (point_val.E.z / E),
                                divE,
                                (dgNode.UP.abs() / U),
                                (dgNode.UP.x / U),
                                (dgNode.UP.y / U),
                                (dgNode.UP.z / U),
                                (dgNode.UP.abs() / vA), // Ma
                                (point_val.UE.abs() / U),
                                (point_val.UE.x / U),
                                (point_val.UE.y / U),
                                (point_val.UE.z / U),
                                (point_val.Te / Te),
                                (c_Wpi));

                        buf.append(line, len);
                    }
                }
    }

    for (const std::string& buf : thread_bufs)
        fwrite(buf.data(), 1, buf.size(), fout);

    fclose(fout);
}

template <class GridT, class DensGridT>
void save_grid_levels(
    const GridT& grid, const DensGridT& dens_grid, PlainFilter& plain_filter, index_t from_level, index_t to_level)
{
    for (index_t lv = from_level; lv != to_level; ++lv) {
        plain_filter.set_level(lv);
        save_subgrid(grid, dens_grid, plain_filter);
    }
}

template <class GridT, class DensGridT>
void save_grid(const GridT& grid, const DensGridT& dens_grid, const std::string& group_name)
{
    print_tm("save_grid: " + group_name);
    using namespace PIC;

    Diagnostics group_diag = save_grid_values;

    if (group_name != "all") {
        const ParticleGroups::ParticleGroup group = ParticleGroups::get_group(group_name);
        group_diag                                = group.diag;
    }

    if (group_diag & save_grid_values) {
        UserGridFilters user_grid_filters;
        for (index_t f = 0; f != user_grid_filters.size(); ++f) {
            UserGridFilter user_grid_filter(user_grid_filters.filter_names().at(f));
            save_subgrid(grid, dens_grid, user_grid_filter);
        }

        const size_t m = 1;

        if (Config::save_grid_x_plains()) {
            PlainFilter grid_filter(group_name, PlainFilter::X, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_x() - m);
        }

        if (Config::save_grid_y_plains()) {
            PlainFilter grid_filter(group_name, PlainFilter::Y, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_y() - m);
        }

        if (Config::save_grid_z_plains()) {
            PlainFilter grid_filter(group_name, PlainFilter::Z, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_z() - m);
        }

        if (Config::save_whole_grid()) {
            SaveAllGrid grid_filter(group_name);
            save_subgrid(grid, dens_grid, grid_filter);
        }
    }
}

template <class NodeT> class GridContainer;
struct Cell;
typedef GridContainer<Cell> Grid;

void save_grid_node(const std::string& prefix, const Grid& grid);
void save_raw_grid(const Grid& grid);
