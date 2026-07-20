#include "particles/particles.h"
#include "particles/particle_groups.h"
#include "particles/marker_particles.h"
#include "core/check_particle.h"
#include "io/io_utilities.h"

#include <fmt/core.h>
#include <algorithm>
#include <ostream>

Particles::Particles() {}

index_t Particles::size() const
{
    return data_.size();
}

index_t Particles::total_size() const
{
    return total_size_;
}

void Particles::use_rank_partition(int rank, int size)
{
    partition_rank_ = rank;
    partition_size_ = size > 1 ? size : 1;
}

void Particles::resize(index_t new_size)
{
    total_size_              = new_size;
    const index_t local_size = partitioned_size(new_size);
    data_.resize(local_size);
    inactive_flags_.resize(local_size, 0);
}

const Particle& Particles::operator[](index_t p) const
{
    return data_[p];
}

Particle& Particles::operator[](index_t p)
{
    return data_.at(p);
}

const Particle& Particles::at(index_t p) const
{
    return data_.at(p);
}

Particle& Particles::at(index_t p)
{
    return data_.at(p);
}

void Particles::set(index_t i, const Particle& part)
{
    if (!accepts_global_index(i))
        return;
    data_.at(local_index(i)) = part;
}

void Particles::set(index_t p, const std::string& group_name, const DblVector& vec_r, const DblVector& vec_v, double ni)
{
    if (!accepts_global_index(p))
        return;

    Particle& particle = data_.at(local_index(p));

    particle.r        = vec_r;
    particle.v        = vec_v;
    const int id      = ParticleGroups::get_id_by_name(group_name);
    particle.group_id = id;
    particle.ni       = ni;
}

void Particles::set(index_t p, int group_id, double x, double y, double z, double vx, double vy, double vz, double ni)
{
    if (!accepts_global_index(p))
        return;

    Particle& particle = data_.at(local_index(p));

    particle.r        = DblVector(x, y, z);
    particle.v        = DblVector(vx, vy, vz);
    particle.group_id = group_id;
    particle.ni       = ni;
}

void Particles::add(const Particle& part)
{
    data_.push_back(part);
    inactive_flags_.push_back(0);
    ++total_size_;
}

void Particles::add(const std::string& group_name, const DblVector& vec_r, const DblVector& vec_v, double ni)
{
    const int id = ParticleGroups::get_id_by_name(group_name);
    data_.push_back(Particle(id, vec_r, vec_v, ni));
    inactive_flags_.push_back(0);
    ++total_size_;
}

void Particles::add(int group_id, double x, double y, double z, double vx, double vy, double vz, double ni)
{
    data_.push_back(Particle(group_id, x, y, z, vx, vy, vz, ni));
    inactive_flags_.push_back(0);
    ++total_size_;
}

void Particles::erase(index_t p)
{
    data_.erase(data_.begin() + p);
    inactive_flags_.erase(inactive_flags_.begin() + p);
    if (total_size_ > 0)
        --total_size_;

    MarkerParticles markers;
    markers.remove_idx(p);
}

void Particles::keep_rank_partition(int rank, int size)
{
    if (size <= 1)
        return;

    ParticlesContainer   kept_particles;
    std::vector<uint8_t> kept_inactive_flags;

    for (index_t p = 0; p < data_.size(); ++p) {
        if ((static_cast<int>(p % static_cast<index_t>(size))) == rank) {
            kept_particles.push_back(data_[p]);
            kept_inactive_flags.push_back(p < inactive_flags_.size() ? inactive_flags_[p] : 0);
        }
    }

    data_.swap(kept_particles);
    inactive_flags_.swap(kept_inactive_flags);
    total_size_ = data_.size();
}

void Particles::remove_later(index_t p)
{
    if (p < inactive_flags_.size())
        inactive_flags_[p] = 1;
}

bool Particles::is_inactive(index_t p) const
{
    if (p >= inactive_flags_.size())
        return false;
    return inactive_flags_[p] != 0;
}

size_t Particles::remove_inactives()
{
    return 0;
}

size_t Particles::inactive_count() const
{
    return static_cast<size_t>(std::count(inactive_flags_.begin(), inactive_flags_.end(), uint8_t {1}));
}

void Particle::to_string(std::string& s) const
{
    s = fmt::format("x = {}\n"
                    "y = {}\n"
                    "z = {}\n"
                    "v_x = {}\n"
                    "v_y = {}\n"
                    "v_z = {}",
                    r.x,
                    r.y,
                    r.z,
                    v.x,
                    v.y,
                    v.z);
}

bool Particles::accepts_global_index(index_t p) const
{
    if (partition_size_ <= 1)
        return true;
    return static_cast<int>(p % static_cast<index_t>(partition_size_)) == partition_rank_;
}

index_t Particles::local_index(index_t p) const
{
    if (partition_size_ <= 1)
        return p;
    return p / static_cast<index_t>(partition_size_);
}

index_t Particles::partitioned_size(index_t total_size) const
{
    if (partition_size_ <= 1)
        return total_size;
    if (partition_rank_ < 0 || partition_rank_ >= partition_size_)
        return 0;
    if (static_cast<index_t>(partition_rank_) >= total_size)
        return 0;
    return (total_size - static_cast<index_t>(partition_rank_) + static_cast<index_t>(partition_size_) - 1)
           / static_cast<index_t>(partition_size_);
}

std::ostream& operator<<(std::ostream& out, const Particle& particle)
{
    return (out << particle.r.x << "\t" << particle.r.y << "\t" << particle.r.z << "\t" << particle.v.x << "\t"
                << particle.v.y << "\t" << particle.v.z << "\t" << particle.ni << std::endl);
}
