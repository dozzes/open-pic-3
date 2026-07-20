#include "lua/call_lua_function.h"
#include "config/config.h"
#include "lua/use_lua.h"
#include "io/io_utilities.h"
#include "particles/particles.h"

#include <fmt/core.h>

bool call_lua_function(const char* func_name)
{
    auto&           lua = get_lua_state();
    sol::state_view sv(static_cast<lua_State*>(lua));

    try {
        sol::protected_function f = sv[func_name];
        if (!f.valid()) {
            fmt::print("Function '{}' not found in \"{}\"\n", func_name, PIC::Config::cfg_script_name());
            return false;
        }

        auto result = f();
        if (!result.valid()) {
            sol::error err = result;
            fmt::print("Lua Runtime Error in '{}': {}\n", func_name, err.what());
            return false;
        }
    } catch (const std::exception& e) {
        fmt::print("C++ Exception: {}\n", e.what());
        return false;
    }
    return true;
}

bool lua_validate_particle(const Particle& particle)
{
    auto& lua = get_lua_state();

    bool ok = true;

    try {
        ok = ScriptBridge::Call<bool>(lua, "validate_particle", particle);
    } catch (const std::exception& e) {
        // Include the underlying reason: "not defined" and "threw an error"
        // need different fixes in the config script.
        print(fmt::format("Function <validate_particle()> failed in \"{}\": {}",
                          PIC::Config::cfg_script_name(), e.what()));
        ok = false;
    }

    return ok;
}
