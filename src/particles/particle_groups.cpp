#include "particles/particle_groups.h"
#include <fmt/core.h>
#include <limits>

using namespace std;

const char* ParticleGroups::all_particles_name = "all";

ParticleGroups::State& ParticleGroups::get_instance()
{
    static State instance;
    return instance;
}

ParticleGroups::ParticleGroup::ParticleGroup(
    const string& name_, int id_, double charge_, double mass_, PIC::Diagnostics diag_)
    : name(name_)
    , id(id_)
    , charge(charge_)
    , mass(mass_)
    , diag(diag_)
{}

bool ParticleGroups::create_group(const string& name, double charge, double mass, PIC::Diagnostics diag)
{
    auto& state = get_instance();

    if (state.groups.find(name) != state.groups.end())
        return false;

    int new_id = state.next_id++;

    auto [it, inserted] = state.groups.emplace(piecewise_construct,
                                               forward_as_tuple(name),
                                               forward_as_tuple(name, new_id, charge, mass, diag));

    if (inserted) {
        if (new_id >= static_cast<int>(state.groups_by_id.size()))
            state.groups_by_id.resize(new_id + 1, nullptr);
        state.groups_by_id[new_id] = &it->second;
    }

    if (mass > state.max_mp)
        state.max_mp = mass;

    return true;
}

void ParticleGroups::set_boundary_kind(const string& name, BoundaryKind boundary_kind)
{
    auto& state = get_instance();
    auto  it    = state.groups.find(name);
    if (it == state.groups.end()) {
        throw invalid_argument(fmt::format("ERROR: Group '{}' not found", name));
    }
    it->second.boundary_kind = boundary_kind;
}

int ParticleGroups::get_id_by_name(const string& name)
{
    if (name == all_particles_name)
        return all_particles_id;

    auto& state = get_instance();
    auto  it    = state.groups.find(name);
    if (it == state.groups.end()) {
        throw invalid_argument(fmt::format("ERROR: Group '{}' not found", name));
    }
    return it->second.id;
}

ParticleGroups::ParticleGroup ParticleGroups::get_group(const string& name)
{
    auto& state = get_instance();
    auto  it    = state.groups.find(name);
    if (it == state.groups.end()) {
        throw invalid_argument(fmt::format("ERROR: Group '{}' not found", name));
    }
    return it->second;
}

ParticleGroups::ParticleGroup ParticleGroups::get_group(int id)
{
    auto& state = get_instance();

    if (id >= 0 && id < static_cast<int>(state.groups_by_id.size()) && state.groups_by_id[id])
        return *state.groups_by_id[id];

    throw invalid_argument(fmt::format("ERROR: ID '{}' not found", id));
}

const std::vector<ParticleGroups::ParticleGroup*>& ParticleGroups::groups_by_id_table()
{
    return get_instance().groups_by_id;
}

double ParticleGroups::max_mp()
{
    return get_instance().max_mp;
}

void ParticleGroups::clear()
{
    auto& state = get_instance();
    state.groups.clear();
    state.groups_by_id.clear(); // vector::clear, not map::clear
    state.next_id = 0;
    state.max_mp  = 0.0;
}

vector<string> ParticleGroups::group_names()
{
    auto&          state = get_instance();
    vector<string> names;
    names.reserve(state.groups.size());
    for (const auto& [name, group] : state.groups)
        names.push_back(name);
    return names;
}

size_t ParticleGroups::get_groups_count()
{
    return get_instance().groups.size();
}

void ParticleGroups::to_stream(ostream& os) const
{
    auto& state = get_instance();

    os << "\n[ OPEN-PIC Particle Groups ]\n";
    os << fmt::format("{:-<55}\n", "");
    os << fmt::format("| {:<12} | {:>3} | {:>10} | {:>10} |\n", "Name", "ID", "Charge", "Mass");
    os << fmt::format("{:-<55}\n", "");

    for (const auto& [name, group] : state.groups) {
        os << fmt::format("| {:<12} | {:>3} | {:>10.4f} | {:>10.4f} |\n", group.name, group.id, group.charge, group.mass);
    }
    os << fmt::format("{:-<55}\n", "");
}

ostream& operator<<(ostream& os, const ParticleGroups& pg)
{
    pg.to_stream(os);
    return os;
}
