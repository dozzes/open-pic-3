#pragma once

#include "grid/grid.h"

namespace PIC {

// Neighbor accessor for one grid cell: names the +/-1 offsets used by the
// Yee-staggered stencils so solver code reads as the discrete operators do.
// Shared by the field solver (curl stencils) and the electron model
// (cell-centred j and div(Ue)).
struct CellStencil
{
    Grid&   grid;
    index_t i;
    index_t j;
    index_t k;

    Cell& c() const { return grid(i, j, k); }

    Cell& xm() const { return grid(i - 1, j, k); }
    Cell& xp() const { return grid(i + 1, j, k); }
    Cell& ym() const { return grid(i, j - 1, k); }
    Cell& yp() const { return grid(i, j + 1, k); }
    Cell& zm() const { return grid(i, j, k - 1); }
    Cell& zp() const { return grid(i, j, k + 1); }

    Cell& xm_yp() const { return grid(i - 1, j + 1, k); }
    Cell& xm_zp() const { return grid(i - 1, j, k + 1); }
    Cell& xp_ym() const { return grid(i + 1, j - 1, k); }
    Cell& xp_zm() const { return grid(i + 1, j, k - 1); }
    Cell& ym_zp() const { return grid(i, j - 1, k + 1); }
    Cell& yp_zm() const { return grid(i, j + 1, k - 1); }

    Cell& xp_yp() const { return grid(i + 1, j + 1, k); }
    Cell& xp_zp() const { return grid(i + 1, j, k + 1); }
    Cell& yp_zp() const { return grid(i, j + 1, k + 1); }
};

} // namespace PIC
