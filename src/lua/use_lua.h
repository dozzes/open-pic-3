#pragma once

#include "lua/inc_lua.h"
#include <sol/sol.hpp>
#include <string>
#include <functional>

sol::state& get_lua_state();

namespace ScriptBridge {
template <typename Ret, typename... Args> Ret Call(sol::state_view lua, const std::string& functionName, Args&&... args)
{
    sol::protected_function func = lua[functionName];

    if (!func.valid()) {
        throw std::runtime_error("Script function not found: " + functionName);
    }

    auto result = func(std::forward<Args>(args)...);

    if (!result.valid()) {
        sol::error err = result;
        throw std::runtime_error(err.what());
    }

    // Use 'if constexpr' (C++17) to handle void return type
    if constexpr (std::is_void_v<Ret>) {
        return;
    } else {
        return result.template get<Ret>();
    }
}
} // namespace ScriptBridge