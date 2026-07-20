#pragma once

#include "core/gather_scatter.h"
#include "grid/grid.h"

namespace PIC::Mpi {

void init(int* argc, char*** argv);
void finalize();
void abort(int error_code);

bool enabled();
int  rank();
int  size();
bool is_root();

void allreduce_grid_density(Grid& grid);
void allreduce_density_grid(DensityGrid& density_grid);

} // namespace PIC::Mpi
