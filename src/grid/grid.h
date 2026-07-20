#pragma once

/**
 * OpenPIC Grid Topology (Staggered C-Grid / Face-Edge centered)
 * -----------------------------------------------------------
 * Coordinate system for a single cell (i, j, k):
 * - NP, Te, eta        --- [Cell Center] (i + 1/2, j + 1/2, k + 1/2)
 *
 * - Ex, UEx, UPx  --- [Face X] (i, j + 1/2, k + 1/2)  <-- Interface between (i) and (i-1)
 *
 * - Ey, UEy, UPy  --- [Face Y] (i + 1/2, j, k + 1/2)  <-- Interface between (j) and (j-1)
 *
 * - Ez, UEz, UPz  --- [Face Z] (i + 1/2, j + 1/2, k)  <-- Interface between (k) and (k-1)
 *
 * - Bx            --- [Edge X] (i + 1/2, j, k)        <-- Center of the edge along X
 *
 * - By            --- [Edge Y] (i, j + 1/2, k)        <-- Center of the edge along Y
 *
 * - Bz            --- [Edge Z] (i, j, k + 1/2)        <-- Center of the edge along Z
 */

#include "config/config.h"
#include "grid/vector_3d.h"
#include "util/aligned_alloc.h"

#include <vector>

#include <stdexcept>
#include <iosfwd>

namespace PIC {

enum CellState { cs_active = 0, cs_absorptive = 1, cs_custom = 2 };

} // namespace PIC

struct Density
{
    Density()
        : UP(0.0, 0.0, 0.0)
        , UP_NP(0.0, 0.0, 0.0)
        , NP(0.0)
    {}
    static constexpr bool has_face_density = true;

    DblVector UP;    // face-centered momentum accumulator: sum of (v * weight) over scattered particles
    DblVector UP_NP; // face-centered weight accumulator: denominator used to normalize UP to a velocity
    double    NP;    // ion number density

    bool operator==(const Density& other) const { return NP == other.NP && UP == other.UP && UP_NP == other.UP_NP; }
};

/************************************************************************/
/* Grid class                                                           */
/************************************************************************/

struct Cell
{
    Cell()
        : NP(0.0)
        , B(0.0, 0.0, 0.0)
        , E(0.0, 0.0, 0.0)
        , UE(0.0, 0.0, 0.0)
        , UP(0.0, 0.0, 0.0)
        , UP_NP(0.0, 0.0, 0.0)
        , Te(0.0)
        , eta(0.0)
        , state_(PIC::cs_active)
    {}

    static constexpr bool has_face_density = true;

    Cell(double           NP_,
         const DblVector& vec_B,
         const DblVector& vec_E,
         const DblVector& vec_UE,
         const DblVector& vec_UP,
         double           Te_,
         double           eta_)
        : NP(NP_)
        , B(vec_B)
        , E(vec_E)
        , UE(vec_UE)
        , UP(vec_UP)
        , UP_NP(0.0, 0.0, 0.0)
        , Te(Te_)
        , eta(eta_)
        , state_(PIC::cs_active)
    {}

    Cell(double NP_,
         double Bx,
         double By,
         double Bz,
         double Ex,
         double Ey,
         double Ez,
         double UEx,
         double UEy,
         double UEz,
         double UPx,
         double UPy,
         double UPz,
         double Te_,
         double eta_)
        : NP(NP_)
        , B(Bx, By, Bz)
        , E(Ex, Ey, Ez)
        , UE(UEx, UEy, UEz)
        , UP(UPx, UPy, UPz)
        , UP_NP(0.0, 0.0, 0.0)
        , Te(Te_)
        , eta(eta_)
        , state_(PIC::cs_active)
    {}

    // Rule of zero: every member is copyable/movable on its own, so the
    // compiler-generated copy/move ctor and assignment (memberwise, including
    // the private state_) are exactly what a hand-written one would do.

    void           set_state(PIC::CellState new_state) { state_ = new_state; }
    PIC::CellState state() const { return state_; }

    bool operator==(const Cell& other) const
    {
        return (this->NP == other.NP && this->B == other.B && this->E == other.E && this->UE == other.UE
                && this->UP == other.UP && this->UP_NP == other.UP_NP && this->Te == other.Te && this->eta == other.eta
                && this->state_ == other.state_);
    }

    double    NP;    // ions density (cm^-3)
    DblVector B;     // magnetic field
    DblVector E;     // electric field
    DblVector UE;    // electron velocity
    DblVector UP;    // ion velocity
    DblVector UP_NP; // transient face-density denominator for UP

    double Te;       // Electron temperature (dynamic field)
    double eta;      // Spitzer resistivity

  private:
    PIC::CellState state_;
};

std::ostream& operator<<(std::ostream& out, const Cell& cell);

// 1. Forward declarations for the template friend system
template <class NodeT> class GridContainer;
template <class NodeT> bool operator==(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2);

template <class NodeT> class GridContainer
{
    // Grant friendship to other template instances and the specific equality operator
    template <class U> friend class GridContainer;
    friend bool operator== <NodeT>(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2);

  public:
    using NodeType = NodeT;

    GridContainer()
        : size_x_(0)
        , size_y_(0)
        , size_z_(0)
        , h_(0.0)
    {}

    GridContainer(index_t sx, index_t sy, index_t sz, double h)
        : size_x_(sx)
        , size_y_(sy)
        , size_z_(sz)
        , h_(h)
    {
        resize(size_x_, size_y_, size_z_);
    }

    // Rule of zero: size_x_/size_y_/size_z_/h_ are plain values and data_ is a
    // std::vector, so the compiler-generated copy/move ctor and assignment
    // are exactly what a hand-written one would do.

    // Zeros density moments (NP, UP, UP_NP) before each push-and-scatter step.
    void reset_current()
    {
#pragma omp parallel for
        for (index_t i = 0; i < data_.size(); ++i) {
            NodeT& node = data_[i];
            node.NP     = 0.0;
            node.UP.x   = 0.0;
            node.UP.y   = 0.0;
            node.UP.z   = 0.0;
            if constexpr (NodeT::has_face_density) {
                node.UP_NP.x = 0.0;
                node.UP_NP.y = 0.0;
                node.UP_NP.z = 0.0;
            }
        }
    }

    template <class DenstyGridType> void add_current(const DenstyGridType& dg)
    {
        using DensityNode = typename DenstyGridType::NodeType;

        for (index_t i = 0; i != data_.size(); ++i) {
            const DensityNode& dg_node = dg.data_[i];
            NodeT&             node    = this->data_[i];

            node.NP += dg_node.NP;
            node.UP += dg_node.UP;
            if constexpr (NodeT::has_face_density && DensityNode::has_face_density) {
                node.UP_NP += dg_node.UP_NP;
            }
        }
    }

    // Reduces per-thread density buffers into this grid in parallel (one OMP thread per cell).
    template <class DenstyGridType> void add_current(const std::vector<DenstyGridType>& local_grids)
    {
        using DensityNode = typename DenstyGridType::NodeType;

        const size_t num_threads = local_grids.size();
        const size_t num_cells   = data_.size();

#pragma omp parallel for schedule(static)
        for (index_t i = 0; i < num_cells; ++i) {
            // Temporary accumulator for thread results
            DensityNode thread_sum;

            for (size_t t = 0; t < num_threads; ++t) {
                const auto& cell = local_grids[t].raw_get(i);
                thread_sum.NP += cell.NP;
                thread_sum.UP += cell.UP;
                if constexpr (DensityNode::has_face_density) {
                    thread_sum.UP_NP += cell.UP_NP;
                }
            }

            // accumulate thread results
            data_[i].NP += thread_sum.NP;
            data_[i].UP += thread_sum.UP;
            if constexpr (NodeT::has_face_density && DensityNode::has_face_density) {
                data_[i].UP_NP += thread_sum.UP_NP;
            }
        }
    }

    // --- Accessors ---
    NodeT& operator()(index_t x, index_t y, index_t z) { return data_[(x * size_y_ + y) * size_z_ + z]; }

    const NodeT& operator()(index_t x, index_t y, index_t z) const { return data_[(x * size_y_ + y) * size_z_ + z]; }

    NodeT& at(index_t x, index_t y, index_t z) { return data_.at((x * size_y_ + y) * size_z_ + z); }

    size_t raw_size() const { return data_.size(); }

    const NodeT& raw_get(size_t i) const { return data_.at(i); }

    NodeT& raw_get(size_t i) { return data_.at(i); }

    const NodeT& at(index_t x, index_t y, index_t z) const { return data_.at((x * size_y_ + y) * size_z_ + z); }

    void set(index_t x, index_t y, index_t z, const NodeT& cell) { at(x, y, z) = cell; }

    void resize(index_t new_size_x, index_t new_size_y, index_t new_size_z)
    {
        size_x_ = new_size_x;
        size_y_ = new_size_y;
        size_z_ = new_size_z;
        data_.resize(size_x_ * size_y_ * size_z_);
    }

    void set_boundary_state(PIC::CellState new_state)
    {
        const index_t nx = size_x();
        const index_t ny = size_y();
        const index_t nz = size_z();

        // Fill planes using a lambda to keep it clean and avoid N^3 traversal
        auto fill_plane = [&](index_t x0, index_t x1, index_t y0, index_t y1, index_t z0, index_t z1) {
            // Explicitly use local copies to satisfy MSVC OpenMP parser
            const index_t end_i = x1;
            const index_t end_j = y1;
            const index_t end_k = z1;

#pragma omp parallel for
            for (index_t i = x0; i <= end_i; ++i) {
                for (index_t j = y0; j <= end_j; ++j) {
                    for (index_t k = z0; k <= end_k; ++k) {
                        (*this)(i, j, k).set_state(new_state);
                    }
                }
            }
        };

        // 1. X-axis boundaries (Left & Right)
        fill_plane(0, 0, 0, ny - 1, 0, nz - 1);
        fill_plane(nx - 1, nx - 1, 0, ny - 1, 0, nz - 1);

        // 2. Y-axis boundaries (Front & Back, avoiding X-corners)
        fill_plane(1, nx - 2, 0, 0, 0, nz - 1);
        fill_plane(1, nx - 2, ny - 1, ny - 1, 0, nz - 1);

        // 3. Z-axis boundaries (Top & Bottom, avoiding X and Y corners)
        fill_plane(1, nx - 2, 1, ny - 2, 0, 0);
        fill_plane(1, nx - 2, 1, ny - 2, nz - 1, nz - 1);
    }

    void   set_step(double h) { h_ = h; }
    double step() const { return h_; }

    double cell_volume() const { return (h_ * h_ * h_); }

    index_t size_x() const { return size_x_; }
    index_t size_y() const { return size_y_; }
    index_t size_z() const { return size_z_; }

    double length_x() const { return ((size_x_ - 1) * h_); }
    double length_y() const { return ((size_y_ - 1) * h_); }
    double length_z() const { return ((size_z_ - 1) * h_); }

  private:
    index_t                                         size_x_;
    index_t                                         size_y_;
    index_t                                         size_z_;
    double                                          h_;
    std::vector<NodeT, AlignedAllocator<NodeT, 64>> data_;
};

// --- External Operator Implementations ---

template <class NodeT> bool operator==(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2)
{
    return (g1.size_x_ == g2.size_x_ && g1.size_y_ == g2.size_y_ && g1.size_z_ == g2.size_z_ && g1.h_ == g2.h_
            && g1.data_ == g2.data_);
}

template <class NodeT> bool operator!=(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2)
{
    return !(g1 == g2);
}

using Grid = GridContainer<Cell>;
