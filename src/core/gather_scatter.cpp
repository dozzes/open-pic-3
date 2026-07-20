#include "core/gather_scatter.h"
#include "config/config.h"
#include "lua/use_lua.h"
#include "io/io_utilities.h"
#include "particles/particle_groups.h"
#include "particles/particles.h"
#include "core/check_particle.h"
#include "lua/call_lua_function.h"

#include <fmt/core.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace PIC {

// Interpolates an edge-centered vector (B) to an arbitrary particle position.
// Each component uses a different (e,f,f)/(f,e,f)/(f,f,e) offset combination
// matching the Yee edge topology; the two weight arrays are precomputed once.
void gather_edge(const Grid& grid, const DblVector& at_point, CellVectorValue val, DblVector& ret_vec)
{
    const double h  = grid.step();
    const double px = at_point.x, py = at_point.y, pz = at_point.z;

    if (Config::scatter_alg() == ScatterAlg::NGP) {
        ret_vec.x = (grid(nearest_grid_index(px, h, EdgeXCentering::x()),
                          nearest_grid_index(py, h, EdgeXCentering::y()),
                          nearest_grid_index(pz, h, EdgeXCentering::z())).*val).x;
        ret_vec.y = (grid(nearest_grid_index(px, h, EdgeYCentering::x()),
                          nearest_grid_index(py, h, EdgeYCentering::y()),
                          nearest_grid_index(pz, h, EdgeYCentering::z())).*val).y;
        ret_vec.z = (grid(nearest_grid_index(px, h, EdgeZCentering::x()),
                          nearest_grid_index(py, h, EdgeZCentering::y()),
                          nearest_grid_index(pz, h, EdgeZCentering::z())).*val).z;
        return;
    }
    if (Config::scatter_alg() == ScatterAlg::TSC) {
        ret_vec.x = gather_vector<EdgeXCentering>(grid, at_point, val, &DblVector::x);
        ret_vec.y = gather_vector<EdgeYCentering>(grid, at_point, val, &DblVector::y);
        ret_vec.z = gather_vector<EdgeZCentering>(grid, at_point, val, &DblVector::z);
        return;
    }

    // Home cell indices for edge (offset=0) and face (offset=0.5) centerings.
    const index_t pi_e = to_idx(px, h, 0.0);
    const index_t pi_f = to_idx(px, h, 0.5);
    const index_t pj_e = to_idx(py, h, 0.0);
    const index_t pj_f = to_idx(py, h, 0.5);
    const index_t pk_e = to_idx(pz, h, 0.0);
    const index_t pk_f = to_idx(pz, h, 0.5);

    // Precompute 1D shape weights: 12 R() calls instead of 72.
    double rx_e[2], rx_f[2], ry_e[2], ry_f[2], rz_e[2], rz_f[2];
    for (int d = 0; d < 2; ++d) {
        rx_e[d] = R(static_cast<double>(pi_e + d) * h - px, h);
        rx_f[d] = R((static_cast<double>(pi_f + d) + 0.5) * h - px, h);
        ry_e[d] = R(static_cast<double>(pj_e + d) * h - py, h);
        ry_f[d] = R((static_cast<double>(pj_f + d) + 0.5) * h - py, h);
        rz_e[d] = R(static_cast<double>(pk_e + d) * h - pz, h);
        rz_f[d] = R((static_cast<double>(pk_f + d) + 0.5) * h - pz, h);
    }

    // Bx: EdgeX centering (x=face/0.5, y=edge/0.0, z=edge/0.0)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_f[di] * ry_e[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_f + di, pj_e + dj, pk_e + dk).*val).x * wxy * rz_e[dk];
            }
        ret_vec.x = sum;
    }

    // By: EdgeY centering (x=edge/0.0, y=face/0.5, z=edge/0.0)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_e[di] * ry_f[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_e + di, pj_f + dj, pk_e + dk).*val).y * wxy * rz_e[dk];
            }
        ret_vec.y = sum;
    }

    // Bz: EdgeZ centering (x=edge/0.0, y=edge/0.0, z=face/0.5)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_e[di] * ry_e[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_e + di, pj_e + dj, pk_f + dk).*val).z * wxy * rz_f[dk];
            }
        ret_vec.z = sum;
    }
}

// Interpolates a face-centered vector (E, UP, UE) to an arbitrary particle position.
// Offset pattern: Ex=(e,f,f), Ey=(f,e,f), Ez=(f,f,e) — dual to gather_edge.
void gather_face(const Grid& grid, const DblVector& at_point, CellVectorValue val, DblVector& ret_vec)
{
    const double h  = grid.step();
    const double px = at_point.x, py = at_point.y, pz = at_point.z;

    if (Config::scatter_alg() == ScatterAlg::NGP) {
        ret_vec.x = (grid(nearest_grid_index(px, h, FaceXCentering::x()),
                          nearest_grid_index(py, h, FaceXCentering::y()),
                          nearest_grid_index(pz, h, FaceXCentering::z())).*val).x;
        ret_vec.y = (grid(nearest_grid_index(px, h, FaceYCentering::x()),
                          nearest_grid_index(py, h, FaceYCentering::y()),
                          nearest_grid_index(pz, h, FaceYCentering::z())).*val).y;
        ret_vec.z = (grid(nearest_grid_index(px, h, FaceZCentering::x()),
                          nearest_grid_index(py, h, FaceZCentering::y()),
                          nearest_grid_index(pz, h, FaceZCentering::z())).*val).z;
        return;
    }
    if (Config::scatter_alg() == ScatterAlg::TSC) {
        ret_vec.x = gather_vector<FaceXCentering>(grid, at_point, val, &DblVector::x);
        ret_vec.y = gather_vector<FaceYCentering>(grid, at_point, val, &DblVector::y);
        ret_vec.z = gather_vector<FaceZCentering>(grid, at_point, val, &DblVector::z);
        return;
    }

    // Home cell indices for edge (offset=0) and face (offset=0.5) centerings.
    const index_t pi_e = to_idx(px, h, 0.0);
    const index_t pi_f = to_idx(px, h, 0.5);
    const index_t pj_e = to_idx(py, h, 0.0);
    const index_t pj_f = to_idx(py, h, 0.5);
    const index_t pk_e = to_idx(pz, h, 0.0);
    const index_t pk_f = to_idx(pz, h, 0.5);

    // Precompute 1D shape weights: 12 R() calls instead of 72.
    double rx_e[2], rx_f[2], ry_e[2], ry_f[2], rz_e[2], rz_f[2];
    for (int d = 0; d < 2; ++d) {
        rx_e[d] = R(static_cast<double>(pi_e + d) * h - px, h);
        rx_f[d] = R((static_cast<double>(pi_f + d) + 0.5) * h - px, h);
        ry_e[d] = R(static_cast<double>(pj_e + d) * h - py, h);
        ry_f[d] = R((static_cast<double>(pj_f + d) + 0.5) * h - py, h);
        rz_e[d] = R(static_cast<double>(pk_e + d) * h - pz, h);
        rz_f[d] = R((static_cast<double>(pk_f + d) + 0.5) * h - pz, h);
    }

    // Ex: FaceX centering (x=edge/0.0, y=face/0.5, z=face/0.5)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_e[di] * ry_f[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_e + di, pj_f + dj, pk_f + dk).*val).x * wxy * rz_f[dk];
            }
        ret_vec.x = sum;
    }

    // Ey: FaceY centering (x=face/0.5, y=edge/0.0, z=face/0.5)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_f[di] * ry_e[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_f + di, pj_e + dj, pk_f + dk).*val).y * wxy * rz_f[dk];
            }
        ret_vec.y = sum;
    }

    // Ez: FaceZ centering (x=face/0.5, y=face/0.5, z=edge/0.0)
    {
        double sum = 0.0;
        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj) {
                const double wxy = rx_f[di] * ry_f[dj];
                for (int dk = 0; dk < 2; ++dk)
                    sum += (grid(pi_f + di, pj_f + dj, pk_e + dk).*val).z * wxy * rz_e[dk];
            }
        ret_vec.z = sum;
    }
}

double gather_center(const Grid& grid, const DblVector& at_point, CellScalarValue val)
{
    return gather_scalar<CellCentering>(grid, at_point, val);
}

void from_grid_to_point(const Grid& grid, const DblVector& vec_to_point, Grid::NodeType& val_at_point)
{
    // interpolate NP
    val_at_point.NP = gather_center(grid, vec_to_point, &Cell::NP);
    val_at_point.Te = gather_center(grid, vec_to_point, &Cell::Te);

    // interpolate B
    gather_edge(grid, vec_to_point, &Cell::B, val_at_point.B);

    // interpolate E
    gather_face(grid, vec_to_point, &Cell::E, val_at_point.E);

    // interpolate UP
    gather_face(grid, vec_to_point, &Cell::UP, val_at_point.UP);

    // interpolate UE
    gather_face(grid, vec_to_point, &Cell::UE, val_at_point.UE);
}

namespace {

// Shared failure contract for every set_boundary_conditions overload: the
// underlying Lua/C++ error is already printed by the time we get here, so the
// exception carries the actionable message (which callback, which script).
// A broken boundary callback must stop the run, not degrade it silently.
[[noreturn]] void throw_boundary_error(const std::string& lua_func_name)
{
    throw std::domain_error(fmt::format("Please check function \"{}\" in \"{}\" and restart.",
                                        lua_func_name, PIC::Config::cfg_script_name()));
}

} // namespace

void set_boundary_conditions(const std::string& lua_func_name)
{
    // call_lua_function prints the underlying error itself.
    if (!call_lua_function(lua_func_name.c_str()))
        throw_boundary_error(lua_func_name);
}

void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg)
{
    try {
        // Pass dg directly, sol2 handles the reference
        ScriptBridge::Call<void>(get_lua_state(), lua_func_name, dg);
    } catch (const std::exception& e) {
        fmt::print("Lua error in '{}': {}\n", lua_func_name, e.what());
        throw_boundary_error(lua_func_name);
    }
}

void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg, const std::string& group_name)
{
    try {
        ScriptBridge::Call<void>(get_lua_state(), lua_func_name, dg, group_name);
    } catch (const std::exception& e) {
        fmt::print("Lua error in '{}': {}\n", lua_func_name, e.what());
        throw_boundary_error(lua_func_name);
    }
}

/**
 * @brief Standard particle-to-grid deposition (Scatter).
 * Matches OpenPIC Staggered C-Grid Topology:
 * - NP (Density)      -> Cell Center (i+0.5, j+0.5, k+0.5)
 * - UP.x (Velocity X) -> Face X      (i,     j+0.5, k+0.5)
 * - UP.y (Velocity Y) -> Face Y      (i+0.5, j,     k+0.5)
 * - UP.z (Velocity Z) -> Face Z      (i+0.5, j+0.5, k)
 */
void scatter_particle_std(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    const double h = Config::h();

    // Pre-calculate total charge and current components.
    // Use the flat id-table to avoid a hash-map lookup + get_instance() on every call.
    const auto&     id_table = ParticleGroups::groups_by_id_table();
    const double    q        = particle.ni * id_table[particle.group_id]->charge;
    const DblVector qv       = {q * particle.v.x, q * particle.v.y, q * particle.v.z};

    // Use new position
    const double px = particle.r.x + dr.x;
    const double py = particle.r.y + dr.y;
    const double pz = particle.r.z + dr.z;

    // Determine the base cell index (left-bottom-near node)
    const index_t pi = to_idx(px, h);
    const index_t pj = to_idx(py, h);
    const index_t pk = to_idx(pz, h);

    // Guard: the 3x3x3 stencil spans [pi-1, pi+1]. If any base index is 0 the
    // lower neighbour would underflow (index_t is unsigned) and access garbage.
    // is_particle_can_scatter() blocks most such particles, but the initial
    // scatter in simulate() calls us directly, so we must check here too.
    if (pi < 1 || pj < 1 || pk < 1 || pi + 1 >= grid.size_x() || pj + 1 >= grid.size_y() || pk + 1 >= grid.size_z())
        return;

    using DensityNode = DensityGrid::NodeType;

    // Precompute 1D shape function values for all stencil positions.
    // Each 3D R(x,y,z) = R(x)*R(y)*R(z), so we precompute the 1D factors
    // (center = half-offset, face = integer-offset) for each axis.
    // 18 R() calls total instead of 108 per particle.
    double rxc[3], rxf[3], ryc[3], ryf[3], rzc[3], rzf[3];
    for (int d = 0; d < 3; ++d) {
        rxc[d] = R((static_cast<double>(pi - 1 + d) + 0.5) * h - px, h);
        rxf[d] = R(static_cast<double>(pi - 1 + d) * h - px, h);
        ryc[d] = R((static_cast<double>(pj - 1 + d) + 0.5) * h - py, h);
        ryf[d] = R(static_cast<double>(pj - 1 + d) * h - py, h);
        rzc[d] = R((static_cast<double>(pk - 1 + d) + 0.5) * h - pz, h);
        rzf[d] = R(static_cast<double>(pk - 1 + d) * h - pz, h);
    }

    // 3x3x3 Cloud-in-Cell (CIC) deposition stencil.
    // Hoist xy cross-products out of the inner k loop for fewer multiplications.
    for (int di = 0; di < 3; ++di)
        for (int dj = 0; dj < 3; ++dj) {
            const double xc_yc = rxc[di] * ryc[dj]; // NP and UP.z share x-center, y-center
            const double xf_yc = rxf[di] * ryc[dj]; // UP.x: x-face, y-center
            const double xc_yf = rxc[di] * ryf[dj]; // UP.y: x-center, y-face
            for (int dk = 0; dk < 3; ++dk) {
                DensityNode& cell = grid(pi - 1 + di, pj - 1 + dj, pk - 1 + dk);
                const double w_np = xc_yc * rzc[dk]; // center, center, center
                const double w_x  = xf_yc * rzc[dk]; // face-x, center, center
                const double w_y  = xc_yf * rzc[dk]; // center, face-y, center
                const double w_z  = xc_yc * rzf[dk]; // center, center, face-z

                cell.NP += q * w_np;
                cell.UP.x += qv.x * w_x;
                cell.UP.y += qv.y * w_y;
                cell.UP.z += qv.z * w_z;
                cell.UP_NP.x += q * w_x;
                cell.UP_NP.y += q * w_y;
                cell.UP_NP.z += q * w_z;
            }
        }
}

/**
 * Nearest-grid-point deposition on the staggered C-grid.  Each scalar or
 * vector component is assigned to the nearest point of its own lattice:
 * NP->[HHH], UPx->[FHH], UPy->[HFH], UPz->[HHF].  Using one common cell for
 * all four quantities would shift the face moments by h/2 and would not be
 * the adjoint of the NGP field gather above.
 */
void scatter_particle_ngp(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    const double h = Config::h();

    const auto&     id_table = ParticleGroups::groups_by_id_table();
    const double    q        = particle.ni * id_table[particle.group_id]->charge;
    const DblVector qv       = {q * particle.v.x, q * particle.v.y, q * particle.v.z};

    const double px = particle.r.x + dr.x;
    const double py = particle.r.y + dr.y;
    const double pz = particle.r.z + dr.z;

    const index_t ci = nearest_grid_index(px, h, CellCentering::x());
    const index_t cj = nearest_grid_index(py, h, CellCentering::y());
    const index_t ck = nearest_grid_index(pz, h, CellCentering::z());

    const index_t xi = nearest_grid_index(px, h, FaceXCentering::x());
    const index_t xj = nearest_grid_index(py, h, FaceXCentering::y());
    const index_t xk = nearest_grid_index(pz, h, FaceXCentering::z());
    const index_t yi = nearest_grid_index(px, h, FaceYCentering::x());
    const index_t yj = nearest_grid_index(py, h, FaceYCentering::y());
    const index_t yk = nearest_grid_index(pz, h, FaceYCentering::z());
    const index_t zi = nearest_grid_index(px, h, FaceZCentering::x());
    const index_t zj = nearest_grid_index(py, h, FaceZCentering::y());
    const index_t zk = nearest_grid_index(pz, h, FaceZCentering::z());

    const auto inside = [&grid](index_t i, index_t j, index_t k) {
        return i < grid.size_x() && j < grid.size_y() && k < grid.size_z();
    };
    if (!inside(ci, cj, ck) || !inside(xi, xj, xk) || !inside(yi, yj, yk) || !inside(zi, zj, zk))
        return;

    grid(ci, cj, ck).NP += q;

    Density& dx = grid(xi, xj, xk);
    dx.UP.x += qv.x;
    dx.UP_NP.x += q;

    Density& dy = grid(yi, yj, yk);
    dy.UP.y += qv.y;
    dy.UP_NP.y += q;

    Density& dz = grid(zi, zj, zk);
    dz.UP.z += qv.z;
    dz.UP_NP.z += q;
}

/**
 * Quadratic triangular-shaped-cloud deposition on the staggered C-grid.
 * Each moment is deposited on the same lattice from which its corresponding
 * field/velocity component is gathered, making the two operations an
 * adjoint pair away from physical boundaries.
 */
void scatter_particle_tsc(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    const double h = Config::h();

    const auto&     id_table = ParticleGroups::groups_by_id_table();
    const double    q        = particle.ni * id_table[particle.group_id]->charge;
    const DblVector qv       = {q * particle.v.x, q * particle.v.y, q * particle.v.z};
    const DblVector position = particle.r + dr;

    const auto deposit_scalar = [&](double ox, double oy, double oz) {
        const auto pi = static_cast<long long>(std::floor(position.x / h - ox + 0.5));
        const auto pj = static_cast<long long>(std::floor(position.y / h - oy + 0.5));
        const auto pk = static_cast<long long>(std::floor(position.z / h - oz + 0.5));
        for (long long i = pi - 1; i <= pi + 1; ++i)
            for (long long j = pj - 1; j <= pj + 1; ++j)
                for (long long k = pk - 1; k <= pk + 1; ++k) {
                    if (i < 0 || j < 0 || k < 0 || i >= static_cast<long long>(grid.size_x())
                        || j >= static_cast<long long>(grid.size_y()) || k >= static_cast<long long>(grid.size_z()))
                        continue;
                    const double w = R_TSC((i + ox) * h - position.x, h)
                                     * R_TSC((j + oy) * h - position.y, h)
                                     * R_TSC((k + oz) * h - position.z, h);
                    grid(static_cast<index_t>(i), static_cast<index_t>(j), static_cast<index_t>(k)).NP += q * w;
                }
    };

    const auto deposit_vector = [&](double ox, double oy, double oz, double DblVector::* comp, double qv_comp) {
        const auto pi = static_cast<long long>(std::floor(position.x / h - ox + 0.5));
        const auto pj = static_cast<long long>(std::floor(position.y / h - oy + 0.5));
        const auto pk = static_cast<long long>(std::floor(position.z / h - oz + 0.5));
        for (long long i = pi - 1; i <= pi + 1; ++i)
            for (long long j = pj - 1; j <= pj + 1; ++j)
                for (long long k = pk - 1; k <= pk + 1; ++k) {
                    if (i < 0 || j < 0 || k < 0 || i >= static_cast<long long>(grid.size_x())
                        || j >= static_cast<long long>(grid.size_y()) || k >= static_cast<long long>(grid.size_z()))
                        continue;
                    const double w = R_TSC((i + ox) * h - position.x, h)
                                     * R_TSC((j + oy) * h - position.y, h)
                                     * R_TSC((k + oz) * h - position.z, h);
                    Density& node = grid(static_cast<index_t>(i), static_cast<index_t>(j), static_cast<index_t>(k));
                    node.UP.*comp += qv_comp * w;
                    node.UP_NP.*comp += q * w;
                }
    };

    deposit_scalar(CellCentering::x(), CellCentering::y(), CellCentering::z());
    deposit_vector(FaceXCentering::x(), FaceXCentering::y(), FaceXCentering::z(), &DblVector::x, qv.x);
    deposit_vector(FaceYCentering::x(), FaceYCentering::y(), FaceYCentering::z(), &DblVector::y, qv.y);
    deposit_vector(FaceZCentering::x(), FaceZCentering::y(), FaceZCentering::z(), &DblVector::z, qv.z);
}

} // end of namespace PIC
