#include "mpi/mpi_runtime.h"

#include <stdexcept>
#include <vector>

#ifdef OPENPIC_ENABLE_MPI
#include <mpi.h>
#endif

namespace PIC::Mpi {
namespace {

bool g_enabled = false;
int  g_rank    = 0;
int  g_size    = 1;

// NP, UP.xyz, UP_NP.xyz — the 7 density fields reduced across ranks each step.
constexpr size_t density_width = 7;

template <class GridT> std::vector<double> pack_density(const GridT& grid)
{
    std::vector<double> buffer(grid.raw_size() * density_width);

    for (size_t i = 0; i < grid.raw_size(); ++i) {
        const auto&  cell = grid.raw_get(i);
        const size_t base = i * density_width;
        buffer[base + 0]  = cell.NP;
        buffer[base + 1]  = cell.UP.x;
        buffer[base + 2]  = cell.UP.y;
        buffer[base + 3]  = cell.UP.z;
        buffer[base + 4]  = cell.UP_NP.x;
        buffer[base + 5]  = cell.UP_NP.y;
        buffer[base + 6]  = cell.UP_NP.z;
    }

    return buffer;
}

template <class GridT> void unpack_density(GridT& grid, const std::vector<double>& buffer)
{
    for (size_t i = 0; i < grid.raw_size(); ++i) {
        auto&        cell = grid.raw_get(i);
        const size_t base = i * density_width;
        cell.NP           = buffer[base + 0];
        cell.UP.x         = buffer[base + 1];
        cell.UP.y         = buffer[base + 2];
        cell.UP.z         = buffer[base + 3];
        cell.UP_NP.x      = buffer[base + 4];
        cell.UP_NP.y      = buffer[base + 5];
        cell.UP_NP.z      = buffer[base + 6];
    }
}

template <class GridT> void allreduce_density(GridT& grid)
{
    if (!g_enabled || g_size == 1)
        return;

#ifdef OPENPIC_ENABLE_MPI
    std::vector<double> send = pack_density(grid);
    std::vector<double> recv(send.size());

    MPI_Allreduce(send.data(), recv.data(), static_cast<int>(send.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    unpack_density(grid, recv);
#else
    throw std::runtime_error("MPI support is not compiled into this OpenPIC build.");
#endif
}

} // namespace

void init(int* argc, char*** argv)
{
#ifdef OPENPIC_ENABLE_MPI
    int provided = 0;
    // FUNNELED: OpenMP threads may call MPI only from the main thread.
    // MULTIPLE is unnecessary here because all MPI calls are in the main thread;
    // requiring it would exclude MS-MPI which only guarantees FUNNELED.
    MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
        throw std::runtime_error("MPI implementation does not provide MPI_THREAD_FUNNELED.");

    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    g_enabled = true;
#else
    (void)argc;
    (void)argv;
    throw std::runtime_error("This OpenPIC build was compiled without MPI support.");
#endif
}

void finalize()
{
#ifdef OPENPIC_ENABLE_MPI
    if (g_enabled) {
        MPI_Finalize();
        g_enabled = false;
        g_rank    = 0;
        g_size    = 1;
    }
#endif
}

void abort(int error_code)
{
#ifdef OPENPIC_ENABLE_MPI
    if (g_enabled)
        MPI_Abort(MPI_COMM_WORLD, error_code);
#endif
}

bool enabled()
{
    return g_enabled;
}

int rank()
{
    return g_rank;
}

int size()
{
    return g_size;
}

bool is_root()
{
    return g_rank == 0;
}

void allreduce_grid_density(Grid& grid)
{
    allreduce_density(grid);
}

void allreduce_density_grid(DensityGrid& density_grid)
{
    allreduce_density(density_grid);
}

} // namespace PIC::Mpi
