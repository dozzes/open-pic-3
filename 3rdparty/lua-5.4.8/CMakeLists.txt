cmake_minimum_required(VERSION 3.10)
project(lua C)

# Ensure CMake looks in the 'src' folder for the files
set(LUA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# 1. CORE SOURCE FILES (Internal logic)
set(CORE_C
    ${LUA_SRC_DIR}/lapi.c ${LUA_SRC_DIR}/lcode.c ${LUA_SRC_DIR}/lctype.c 
    ${LUA_SRC_DIR}/ldebug.c ${LUA_SRC_DIR}/ldo.c ${LUA_SRC_DIR}/ldump.c 
    ${LUA_SRC_DIR}/lfunc.c ${LUA_SRC_DIR}/lgc.c ${LUA_SRC_DIR}/llex.c 
    ${LUA_SRC_DIR}/lmem.c ${LUA_SRC_DIR}/lobject.c ${LUA_SRC_DIR}/lopcodes.c 
    ${LUA_SRC_DIR}/lparser.c ${LUA_SRC_DIR}/lstate.c ${LUA_SRC_DIR}/lstring.c 
    ${LUA_SRC_DIR}/ltable.c ${LUA_SRC_DIR}/ltm.c ${LUA_SRC_DIR}/lundump.c 
    ${LUA_SRC_DIR}/lvm.c ${LUA_SRC_DIR}/lzio.c
)

# 2. LIBRARY FILES (Standard Libraries)
set(LIB_C
    ${LUA_SRC_DIR}/lauxlib.c ${LUA_SRC_DIR}/lbaselib.c ${LUA_SRC_DIR}/lcorolib.c 
    ${LUA_SRC_DIR}/ldblib.c ${LUA_SRC_DIR}/liolib.c ${LUA_SRC_DIR}/lmathlib.c 
    ${LUA_SRC_DIR}/loadlib.c ${LUA_SRC_DIR}/loslib.c ${LUA_SRC_DIR}/lstrlib.c 
    ${LUA_SRC_DIR}/ltablib.c ${LUA_SRC_DIR}/lutf8lib.c ${LUA_SRC_DIR}/linit.c
)

# 3. CREATE THE STATIC LIBRARY
# Note: We exclude lua.c and luac.c here so they don't clash with your main app
add_library(liblua STATIC ${CORE_C} ${LIB_C})

# Define the include directory so sol2 can find lua.h
target_include_directories(liblua PUBLIC ${LUA_SRC_DIR})

# 4. OPTIONAL: CREATE STANDALONE EXECUTABLES
# Only build these if you want the lua.exe interpreter
add_executable(lua_shell ${LUA_SRC_DIR}/lua.c)
target_link_libraries(lua_shell PRIVATE liblua)

add_executable(lua_compiler ${LUA_SRC_DIR}/luac.c)
target_link_libraries(lua_compiler PRIVATE liblua)