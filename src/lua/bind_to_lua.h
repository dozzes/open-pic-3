#pragma once

#include "util/opic_fwd.h"

bool bind_to_lua(sol::state& lua, const char* lua_cfg_file_name, Grid& grid, Particles& particles);
