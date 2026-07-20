#include "grid/grid_filters.h"
#include "config/config.h"
#include "lua/use_lua.h"
#include "io/io_utilities.h"

#include <fmt/core.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

// static
vector<string> UserGridFilters::filters_ = vector<string>();

UserGridFilter::UserGridFilter(const string& name)
    : GridFilter(name)
{}

bool UserGridFilter::operator()(const DblVector& node, const DblVector& pos) const
{
    auto& lua = get_lua_state();

    bool do_save = false;

    try {
        do_save = ScriptBridge::Call<bool>(lua, name_, std::ref(node), std::ref(pos));
    } catch (exception& e) {
        cerr << "\n" << name_ << " : " << e.what() << endl;
    }

    return do_save;
}

SaveAllGrid::SaveAllGrid(const string& grid_group_name)
    : GridFilter(grid_group_name)
{
    name_ += "_xyz";
}

const string PlainFilter::plain_tags_ = "xyz";

PlainFilter::PlainFilter(const string& grid_group_name, Plain plain, index_t level)
    : GridFilter("")
    , plain_(plain)
    , level_(level)
{
    name_ = fmt::format("{}_{}_{}", grid_group_name, plain_tags_[plain_], level_);
}

void PlainFilter::set_level(index_t level)
{
    level_     = level;
    size_t pos = name_.find_last_of('_');
    name_      = fmt::format("{}{}", name_.substr(0, pos + 1), level_);
}

bool PlainFilter::operator()(const DblVector& node, const DblVector& /*pos*/) const
{
    return (level_ == (plain_ == X ? node.x : plain_ == Y ? node.y : node.z));
}
