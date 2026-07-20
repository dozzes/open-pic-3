#include "lua/use_lua.h"

// Usage:
//      auto& lua = get_lua_state();
sol::state& get_lua_state()
{
    static sol::state lua;
    return lua;
}
