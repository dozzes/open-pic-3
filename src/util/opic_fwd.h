#pragma once

#include <iosfwd>
#include <cmath>
#include <vector>

template <class NodeT> class GridContainer;

struct Cell;
typedef GridContainer<Cell> Grid;

template <class ScalarT> class Vector3D;

typedef Vector3D<double> DblVector;

typedef double(DblVector::* VectorComp);

struct Particle;

class Particles;

typedef size_t index_t;

inline index_t to_idx(double pos, double h, double offset = 0.0)
{
    return static_cast<index_t>(std::floor(pos / h - offset));
}

namespace sol {

class state;
};
