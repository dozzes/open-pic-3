#pragma once

#include "grid/vector_3d.h"
#include "config/config.h"
#include "config/constants.h"

#include <omp.h>
#include <cmath>
#include <random>

// Per-thread deterministic engine, seeded lazily on first use from
// pic_parameters.random_seed. The OpenMP thread id is mixed in (golden-ratio
// stride) so concurrent threads draw independent streams; thread 0 uses the
// seed exactly. Because seeding happens at first use, changing random_seed
// after any rnd() call has no effect on that thread's stream.
//
// MPI: Lua-side particle placement (sim/lib/init_particles.lua) runs on every
// rank with the same seed and consumes the stream identically, so all ranks
// generate the same global particle set; the rank partition then just selects
// disjoint subsets. Any future C++ placement loop must stay single-threaded
// for the same reason -- a parallel loop maps particle index to RNG stream
// via thread id, which is chunking-dependent and breaks this invariant.
inline std::mt19937_64& rnd_engine()
{
    thread_local std::mt19937_64 engine(
        PIC::Config::random_seed()
        + 0x9E3779B97F4A7C15ULL * static_cast<unsigned long long>(omp_get_thread_num()));
    return engine;
}

inline double rnd(double r_min = 0.0, double r_max = 1.0)
{
    return std::uniform_real_distribution<double>(r_min, r_max)(rnd_engine());
}

// Uniform point in a ball of radius R: r ~ R*u^(1/3) makes the radial CDF
// proportional to volume, cos(theta) uniform in [-1,1] removes pole clustering.
inline void rnd_ball(double R, DblVector& result)
{
    double r = R * pow(rnd(), 1.0 / 3.0);
    double t = acos(2.0 * rnd() - 1);
    double f = 2.0 * PIC::Constants::pi() * rnd();

    result.x = r * sin(t) * cos(f);
    result.y = r * sin(t) * sin(f);
    result.z = r * cos(t);
}
