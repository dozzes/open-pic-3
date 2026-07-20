/*
Full-integer coordinates [F] are at i*h, j*h, k*h.
Half-integer coordinates [H] are at (i + 1/2)*h, (j + 1/2)*h, (k + 1/2)*h.

Charge density NP is defined at [HHH].
Ex, Ey, Ez and Jx, Jy, Jz are defined at [FHH], [HFH], [HHF].
Bx, By, Bz are defined at [HFF], [FHF], [FFH].
*/

#pragma once

#include "grid/grid.h"
#include "particles/particles.h"
#include "config/config.h"

#include <cmath>
#include <vector>
#include <string>

namespace PIC {

/***************************************************************************
 * 1D particle form-factor                                                  *
 ***************************************************************************/
inline double R(double x, double h)
{
    return ((fabs(x) < h) ? (1.0 - fabs(x) / h) : 0.0);
}

// Quadratic B-spline (triangular-shaped-cloud, TSC) particle form factor.
// Its support spans the nearest lattice point and its two neighbours.
inline double R_TSC(double x, double h)
{
    const double u = std::fabs(x / h);
    if (u < 0.5)
        return 0.75 - u * u;
    if (u < 1.5) {
        const double d = 1.5 - u;
        return 0.5 * d * d;
    }
    return 0.0;
}

/***************************************************************************
 * 3D particle form-factor                                                  *
 ***************************************************************************/
inline double R(double x, double y, double z, double h)
{
    return R(x, h) * R(y, h) * R(z, h);
}

/***************************************************************************
 * Pointer to Cell value                                                    *
 ***************************************************************************/
using CellVectorValue = DblVector Cell::*;
using CellScalarValue = double Cell::*;

// Nearest point on a lattice whose samples are located at (i+offset)*h.
// Ties are assigned to the upper index, matching the usual floor(x+1/2)
// NGP convention. Callers only use this for particles already inside the
// guarded simulation domain, so the unsigned index cannot underflow.
inline index_t nearest_grid_index(double x, double h, double offset)
{
    return static_cast<index_t>(std::floor(x / h - offset + 0.5));
}

/*****************************************************************************
 * Gather edge-centered values.                                               *
 *    at_point - specified position                                           *
 *    val_x, val_y, val_z - Cell required values                              *
 *    ret_vec - result gathered value                                         *
 * Used to interpolate Bx, By, Bz for specified position                      *
 *****************************************************************************/
void gather_edge(const Grid& grid, const DblVector& at_point, CellVectorValue val, DblVector& ret_vec);

/*****************************************************************************
 * Gather face-centered values: Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz      *
 *    at_point - specified position                                           *
 *    val_x, val_y, val_z - Pointers to cell members required                 *
 *    ret_vec - result gathered value                                         *
 * Used to interpolate Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz values        *
 * for specified position.                                                    *
 *****************************************************************************/
void gather_face(const Grid& grid, const DblVector& at_point, CellVectorValue val, DblVector& ret_vec);

/*****************************************************************************
 * Gather cell-centered value.                                                *
 *   at_point - specified position                                            *
 *   val - Pointer to cell member required                                    *
 *   ret_val - result gathered value                                          *
 * Used to interpolate density NP values for specified position.              *
 *                                                                            *
 *****************************************************************************/
double gather_center(const Grid& grid, const DblVector& at_point, CellScalarValue val);

/***************************************************************************
 * interpolates grid values to specified point <vec_to_point>               *
 ***************************************************************************/
void from_grid_to_point(const Grid& grid, const DblVector& vec_to_point, Grid::NodeType& val_at_point);

struct FaceXCentering
{
    static constexpr double x() { return 0.0; }
    static constexpr double y() { return 0.5; }
    static constexpr double z() { return 0.5; }
};

struct FaceYCentering
{
    static constexpr double x() { return 0.5; }
    static constexpr double y() { return 0.0; }
    static constexpr double z() { return 0.5; }
};

struct FaceZCentering
{
    static constexpr double x() { return 0.5; }
    static constexpr double y() { return 0.5; }
    static constexpr double z() { return 0.0; }
};

struct EdgeXCentering
{
    static constexpr double x() { return 0.5; }
    static constexpr double y() { return 0.0; }
    static constexpr double z() { return 0.0; }
};

struct EdgeYCentering
{
    static constexpr double x() { return 0.0; }
    static constexpr double y() { return 0.5; }
    static constexpr double z() { return 0.0; }
};

struct EdgeZCentering
{
    static constexpr double x() { return 0.0; }
    static constexpr double y() { return 0.0; }
    static constexpr double z() { return 0.5; }
};

struct CellCentering
{
    static constexpr double x() { return 0.5; }
    static constexpr double y() { return 0.5; }
    static constexpr double z() { return 0.5; }
};

struct NodeCentering
{
    static constexpr double x() { return 0.0; }
    static constexpr double y() { return 0.0; }
    static constexpr double z() { return 0.0; }
};

template <class Centering, class GridT>
double gather_vector(const GridT& grid, const DblVector& at_point, DblVector GridT::NodeType::* vec, double DblVector::* comp)
{
    const double h = grid.step();

    if (Config::scatter_alg() == ScatterAlg::NGP) {
        const index_t i = nearest_grid_index(at_point.x, h, Centering::x());
        const index_t j = nearest_grid_index(at_point.y, h, Centering::y());
        const index_t k = nearest_grid_index(at_point.z, h, Centering::z());
        return grid(i, j, k).*vec.*comp;
    }

    if (Config::scatter_alg() == ScatterAlg::TSC) {
        const auto pi = static_cast<long long>(std::floor(at_point.x / h - Centering::x() + 0.5));
        const auto pj = static_cast<long long>(std::floor(at_point.y / h - Centering::y() + 0.5));
        const auto pk = static_cast<long long>(std::floor(at_point.z / h - Centering::z() + 0.5));
        double ret_val = 0.0;
        for (long long i = pi - 1; i <= pi + 1; ++i)
            for (long long j = pj - 1; j <= pj + 1; ++j)
                for (long long k = pk - 1; k <= pk + 1; ++k) {
                    if (i < 0 || j < 0 || k < 0 || i >= static_cast<long long>(grid.size_x())
                        || j >= static_cast<long long>(grid.size_y()) || k >= static_cast<long long>(grid.size_z()))
                        continue;
                    ret_val += (grid(static_cast<index_t>(i), static_cast<index_t>(j), static_cast<index_t>(k)).*vec.*comp)
                               * R_TSC((i + Centering::x()) * h - at_point.x, h)
                               * R_TSC((j + Centering::y()) * h - at_point.y, h)
                               * R_TSC((k + Centering::z()) * h - at_point.z, h);
                }
        return ret_val;
    }

    // home cell node indexes
    const index_t pi = to_idx(at_point.x, h, Centering::x());
    const index_t pj = to_idx(at_point.y, h, Centering::y());
    const index_t pk = to_idx(at_point.z, h, Centering::z());

    double ret_val = 0.0;

    for (index_t i = pi; i != pi + 2; ++i)
        for (index_t j = pj; j != pj + 2; ++j)
            for (index_t k = pk; k != pk + 2; ++k) {
                ret_val += (grid(i, j, k).*vec.*comp)
                           * R((i + Centering::x()) * h - at_point.x,
                               (j + Centering::y()) * h - at_point.y,
                               (k + Centering::z()) * h - at_point.z,
                               h);
            }

    return ret_val;
}

template <class Centering, class GridT>
double gather_vector(
    const GridT& grid, double x, double y, double z, DblVector GridT::NodeType::* vec, double DblVector::* comp)
{
    return gather_vector<Centering>(grid, DblVector(x, y, z), vec, comp);
}

template <class GridT>
void gather_face(const GridT& grid, const DblVector& at_point, DblVector GridT::NodeType::* val, DblVector& ret_vec)
{
    ret_vec.x = gather_vector<PIC::FaceXCentering>(grid, at_point, val, &DblVector::x);
    ret_vec.y = gather_vector<PIC::FaceYCentering>(grid, at_point, val, &DblVector::y);
    ret_vec.z = gather_vector<PIC::FaceZCentering>(grid, at_point, val, &DblVector::z);
}

template <class Centering, class GridT>
double gather_scalar(const GridT& grid, const DblVector& at_point, double GridT::NodeType::* val)
{
    const double h = grid.step();

    if (Config::scatter_alg() == ScatterAlg::NGP) {
        const index_t i = nearest_grid_index(at_point.x, h, Centering::x());
        const index_t j = nearest_grid_index(at_point.y, h, Centering::y());
        const index_t k = nearest_grid_index(at_point.z, h, Centering::z());
        return grid(i, j, k).*val;
    }


    if (Config::scatter_alg() == ScatterAlg::TSC) {
        const auto pi = static_cast<long long>(std::floor(at_point.x / h - Centering::x() + 0.5));
        const auto pj = static_cast<long long>(std::floor(at_point.y / h - Centering::y() + 0.5));
        const auto pk = static_cast<long long>(std::floor(at_point.z / h - Centering::z() + 0.5));
        double ret_val = 0.0;
        for (long long i = pi - 1; i <= pi + 1; ++i)
            for (long long j = pj - 1; j <= pj + 1; ++j)
                for (long long k = pk - 1; k <= pk + 1; ++k) {
                    if (i < 0 || j < 0 || k < 0 || i >= static_cast<long long>(grid.size_x())
                        || j >= static_cast<long long>(grid.size_y()) || k >= static_cast<long long>(grid.size_z()))
                        continue;
                    ret_val += grid(static_cast<index_t>(i), static_cast<index_t>(j), static_cast<index_t>(k)).*val
                               * R_TSC((i + Centering::x()) * h - at_point.x, h)
                               * R_TSC((j + Centering::y()) * h - at_point.y, h)
                               * R_TSC((k + Centering::z()) * h - at_point.z, h);
                }
        return ret_val;
    }

    // home cell node indexes
    const index_t pi = to_idx(at_point.x, h, Centering::x());
    const index_t pj = to_idx(at_point.y, h, Centering::y());
    const index_t pk = to_idx(at_point.z, h, Centering::z());

    double ret_val = 0.0;

    for (index_t i = pi; i != pi + 2; ++i)
        for (index_t j = pj; j != pj + 2; ++j)
            for (index_t k = pk; k != pk + 2; ++k) {
                ret_val += (grid(i, j, k).*val)
                           * R((i + Centering::x()) * h - at_point.x,
                               (j + Centering::y()) * h - at_point.y,
                               (k + Centering::z()) * h - at_point.z,
                               h);
            }

    return ret_val;
}

using DensityGrid = GridContainer<Density>;

void scatter_particle_std(const Particle& particlep, DensityGrid& grid, const DblVector& dr);
void scatter_particle_ngp(const Particle& particle, DensityGrid& grid, const DblVector& dr);
void scatter_particle_tsc(const Particle& particle, DensityGrid& grid, const DblVector& dr);

/***************************************************************************
 * Set boundary conditions                                                  *
 ***************************************************************************/
void set_boundary_conditions(const std::string& lua_func_name);
void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg);
void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg, const std::string& group_name);

} // namespace PIC
