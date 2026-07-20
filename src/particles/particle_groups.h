#pragma once

#include "util/opic_fwd.h"
#include "config/diagnostics.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>

class ParticleGroups
{
  public:
    static const char* all_particles_name;
    static const int   all_particles_id = -1;

    enum class BoundaryKind { Default, FlowBackground };

    struct ParticleGroup
    {
        ParticleGroup(const std::string& name_, int id_, double charge_, double mass_, PIC::Diagnostics diag_);

        std::string      name;
        int              id;
        double           charge;
        double           mass;
        PIC::Diagnostics diag;
        BoundaryKind     boundary_kind = BoundaryKind::Default;
    };

    ParticleGroups() = default;

    bool create_group(const std::string& name, double charge, double mass, PIC::Diagnostics diag);
    void set_boundary_kind(const std::string& name, BoundaryKind boundary_kind);
    void clear();

    static size_t        get_groups_count();
    static int           get_id_by_name(const std::string& name);
    static ParticleGroup get_group(const std::string& name);
    static ParticleGroup get_group(int id);
    // Returns the flat id→group* table for O(1) access without per-call get_instance() overhead.
    // Cache this pointer before a parallel loop and index directly by group_id.
    static const std::vector<ParticleGroup*>& groups_by_id_table();

    static double                   max_mp();
    static std::vector<std::string> group_names();

    void to_stream(std::ostream& os) const;

  private:
    // Using unordered_map for O(1) string lookup
    typedef std::unordered_map<std::string, ParticleGroup> MapNameToGroup;

    struct State
    {
        MapNameToGroup              groups;
        std::vector<ParticleGroup*> groups_by_id; // indexed by id — O(1) vs hash map
        int                         next_id = 0;
        double                      max_mp  = 0.0;
    };

    static State& get_instance();
};

std::ostream& operator<<(std::ostream& os, const ParticleGroups& pg);
