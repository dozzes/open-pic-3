#pragma once

// Macroparticle container and the Particles collection.
// Units follow the CGS-Gaussian convention used throughout OpenPIC.

#include "util/opic_fwd.h"
#include "grid/vector_3d.h"
#include "util/aligned_alloc.h"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <deque>
#include <vector>

struct alignas(64) Particle
{
    Particle()
        : r(0.0, 0.0, 0.0)
        , v(0.0, 0.0, 0.0)
        , ni(0.0)
        , group_id(-2)
    {}

    Particle(int group_id_, double x_, double y_, double z_, double v_x_, double v_y_, double v_z_, double ni_)
        : r(x_, y_, z_)
        , v(v_x_, v_y_, v_z_)
        , ni(ni_)
        , group_id(group_id_)
    {}

    Particle(int group_id_, const DblVector& vec_r, const DblVector& vec_v, double ni_)
        : r(vec_r)
        , v(vec_v)
        , ni(ni_)
        , group_id(group_id_)
    {}

    void to_string(std::string& s) const;

    DblVector r;
    DblVector v;
    double    ni;       // number of real plasma ions this macroparticle represents

    int group_id;
};

std::ostream& operator<<(std::ostream& out, const Particle& p);

class Particles
{
    typedef std::deque<Particle, AlignedAllocator<Particle, 64>> ParticlesContainer;

  public:
    Particles();

    // size()       — particles held by this rank (after partition).
    // total_size() — total across all ranks, as recorded by the Lua script.
    index_t size() const;
    index_t total_size() const;
    void    resize(index_t new_size);

    // Must be called before bind_to_lua so injection callbacks see the correct
    // process_idx and skip particle creation on non-owner ranks.
    void use_rank_partition(int rank, int size);

    const Particle& operator[](index_t p) const;
    Particle&       operator[](index_t p);

    const Particle& at(index_t p) const;
    Particle&       at(index_t p);

    // Deferred removal: remove_later() flags without invalidating iterators.
    // remove_inactives() intentionally does not compact now: compaction changes
    // particle indices and breaks the current OpenMP particle loops.
    size_t remove_inactives();
    size_t inactive_count() const;
    void   remove_later(index_t p);
    bool   is_inactive(index_t p) const;

    // for using from Lua
    void set(index_t p, const Particle& part);

    void set(index_t p, const std::string& group_name, const DblVector& vec_r, const DblVector& vec_v, double ni);

    // Fast path for bulk init from Lua: caller resolves group_id once and
    // passes raw components, skipping the per-call group-name lookup and
    // DblVector userdata construction that set(..., group_name, ...) incurs.
    void set(index_t p, int group_id, double x, double y, double z, double vx, double vy, double vz, double ni);

    void add(const Particle& part);

    void add(const std::string& group_name, const DblVector& vec_r, const DblVector& vec_v, double ni);

    // Fast path counterpart to the raw-double set() overload above.
    void add(int group_id, double x, double y, double z, double vx, double vy, double vz, double ni);

    void erase(index_t p);

    // Called after bind_to_lua to drop particles not owned by this rank.
    // Ownership rule: particle p belongs to rank (p % size).
    void keep_rank_partition(int rank, int size);

    ParticlesContainer::const_iterator begin() const { return data_.begin(); }
    ParticlesContainer::const_iterator end() const { return data_.end(); }

    ParticlesContainer::iterator begin() { return data_.begin(); }
    ParticlesContainer::iterator end() { return data_.end(); }

  private:
    bool    accepts_global_index(index_t p) const;
    index_t local_index(index_t p) const;
    index_t partitioned_size(index_t total_size) const;

    ParticlesContainer   data_;
    std::vector<uint8_t> inactive_flags_; // one flag per particle, 1 = marked for removal
    index_t              total_size_     = 0;
    int                  partition_rank_ = 0;
    int                  partition_size_ = 1;
};
