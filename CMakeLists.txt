cmake_minimum_required(VERSION 3.13)
project(open-pic LANGUAGES CXX C)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "" FORCE)

# 1. GLOBAL SETTINGS
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Mute "unsafe" function warnings in MSVC
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# 2. ADD SUBDIRECTORIES
add_subdirectory(3rdparty/lua-5.4.8)

option(OPENPIC_BUILD_TESTS "Build OpenPIC verification tests" OFF)
option(OPENPIC_PROFILE_PARTICLE_PHASES "Print detailed move_particles_half_time phase timings" OFF)
option(OPENPIC_PARTICLE_STATIC_SCHEDULE "Use OpenMP static schedule in move_particles_half_time" OFF)

# 3. DEFINE THE EXECUTABLE
# Explicit source list (no GLOB): a new file that is not registered here must
# fail the build loudly instead of being silently compiled against a stale
# CMake cache ("собралось старое").
set(SOURCES
    src/config/config.cpp
    src/core/check_particle.cpp
    src/core/electron_model.cpp
    src/core/field_solver.cpp
    src/core/filters.cpp
    src/core/gather_scatter.cpp
    src/core/sim_loop.cpp
    src/grid/grid.cpp
    src/grid/grid_filters.cpp
    src/io/io_utilities.cpp
    src/io/save_grid.cpp
    src/io/save_particles.cpp
    src/lua/bind_to_lua.cpp
    src/lua/call_lua_function.cpp
    src/lua/use_lua.cpp
    src/mpi/mpi_runtime.cpp
    src/particles/marker_particles.cpp
    src/particles/particle_groups.cpp
    src/particles/particles.cpp
)
set(HEADERS
    src/config/config.h
    src/config/constants.h
    src/config/diagnostics.h
    src/core/cell_stencil.h
    src/core/check_particle.h
    src/core/electron_model.h
    src/core/field_solver.h
    src/core/filters.h
    src/core/gather_scatter.h
    src/core/heat_solver.h
    src/core/simulate.h
    src/core/solver_workspace.h
    src/grid/grid.h
    src/grid/grid_filters.h
    src/grid/vector_3d.h
    src/grid/vector_3d.inl
    src/io/io_utilities.h
    src/io/save_grid.h
    src/io/save_particles.h
    src/lua/bind_to_lua.h
    src/lua/call_lua_function.h
    src/lua/inc_lua.h
    src/lua/use_lua.h
    src/mpi/mpi_runtime.h
    src/particles/marker_particles.h
    src/particles/particle_groups.h
    src/particles/particles.h
    src/util/aligned_alloc.h
    src/util/msvc.h
    src/util/opic_fwd.h
    src/util/utility.h
)

add_library(openpic_core OBJECT ${SOURCES} ${HEADERS})
add_executable(${PROJECT_NAME} "src/main.cpp" $<TARGET_OBJECTS:openpic_core>)
source_group(TREE "${CMAKE_SOURCE_DIR}/src" FILES ${SOURCES} ${HEADERS})

# 4. PATHS
set(SOL2_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/3rdparty/sol2-3.3.0/include")

# 5. OPTIMIZATION AND COMPILER FLAGS
# Single carrier for flags/includes/definitions: every target that compiles
# project headers must build under identical semantics, or tests validate
# different codegen than production (this drift already happened once).
add_library(openpic_build_flags INTERFACE)

if(MSVC)
    target_compile_options(openpic_build_flags INTERFACE
        /W3 /GR /EHsc /permissive- /Zc:__cplusplus /Zc:twoPhase- /arch:AVX2 /bigobj
    )

    # Debug: no optimizations, full checks
    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:Debug>:/Od;/Zi;/RTC1;/D_DEBUG>"
    )

    # Release: maximum performance
    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:Release>:/Ox;/GL;/fp:precise;/DNDEBUG>"
    )

    # RelWithDebInfo: optimized but debuggable (NO /GL!)
    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:RelWithDebInfo>:/Ox;/Oy-;/fp:precise;/DNDEBUG;/Zi>"
    )

    # /GL objects require /LTCG on every executable that links them.
    target_link_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:Release>:/LTCG>"
        "$<$<CONFIG:RelWithDebInfo>:/DEBUG>"
    )

else()
    # GCC/Clang
    target_compile_options(openpic_build_flags INTERFACE -Wall -Wextra)

    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:Debug>:-Og;-g>"
    )

    # -march=native maximizes throughput but makes results machine-dependent;
    # V&V runs that must reproduce across hosts disable it.
    option(OPENPIC_NATIVE_ARCH "Compile with -march=native (machine-dependent codegen)" ON)
    if(OPENPIC_NATIVE_ARCH)
        set(OPENPIC_ARCH_FLAG "-march=native")
    else()
        set(OPENPIC_ARCH_FLAG "")
    endif()

    # Release: max perf
    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:Release>:-O3;${OPENPIC_ARCH_FLAG};-DNDEBUG>"
    )

    # RelWithDebInfo: optimized + debuggable (NO fast-math!)
    target_compile_options(openpic_build_flags INTERFACE
        "$<$<CONFIG:RelWithDebInfo>:-O3;${OPENPIC_ARCH_FLAG};-g;-DNDEBUG>"
    )
endif()

# 6. INCLUDE DIRECTORIES
target_include_directories(openpic_build_flags INTERFACE
    "${CMAKE_SOURCE_DIR}/src"
    "${SOL2_INCLUDE_DIR}"
    "${CMAKE_SOURCE_DIR}/3rdparty/lua-5.4.8/src"
    "${CMAKE_SOURCE_DIR}/3rdparty/pocketfft"
    "${CMAKE_SOURCE_DIR}/3rdparty/fmt-10.2.1/include"
)
target_compile_definitions(openpic_build_flags INTERFACE FMT_HEADER_ONLY)

target_link_libraries(openpic_core PRIVATE openpic_build_flags)

if(OPENPIC_PROFILE_PARTICLE_PHASES)
    target_compile_definitions(openpic_core PRIVATE OPENPIC_PROFILE_PARTICLE_PHASES)
endif()

if(OPENPIC_PARTICLE_STATIC_SCHEDULE)
    target_compile_definitions(openpic_core PRIVATE OPENPIC_PARTICLE_STATIC_SCHEDULE)
endif()

# 7. LINK LIBRARIES
target_link_libraries(${PROJECT_NAME} PRIVATE
    liblua
    openpic_build_flags
)

# 8. OPENMP SUPPORT (Using LLVM backend for MSVC 2022)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    if(MSVC)
        # Use LLVM-based OpenMP for better performance and C++11/14/17 support
        target_compile_options(openpic_build_flags INTERFACE /openmp:llvm)
        # Explicitly link the OpenMP runtime for LLVM-MSVC
        target_link_libraries(openpic_build_flags INTERFACE libomp)
    else()
        target_link_libraries(openpic_build_flags INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

# 8b. OPTIONAL MPI SUPPORT
option(OPENPIC_ENABLE_MPI "Enable MPI distributed execution mode" OFF)
if(OPENPIC_ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(openpic_build_flags INTERFACE MPI::MPI_CXX)
    target_compile_definitions(openpic_build_flags INTERFACE OPENPIC_ENABLE_MPI)
endif()

if(OPENPIC_BUILD_TESTS)
    enable_testing()

    add_executable(openpic_tests "tests/openpic_tests.cpp" $<TARGET_OBJECTS:openpic_core>)
    target_link_libraries(openpic_tests PRIVATE liblua openpic_build_flags)

    add_test(NAME openpic_tests COMMAND openpic_tests)
endif()

# 9. OUTPUT AND DEBUGGER PROPERTIES
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
if(OPENPIC_BUILD_TESTS)
    # Keep both executables in the same bin/ directory so run instructions
    # and scripts need only one path convention.
    set_target_properties(openpic_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
endif()

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES
        VS_DEBUGGER_COMMAND_ARGUMENTS "main.lua"
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        # KMP_BLOCKTIME=0: OpenMP worker threads yield immediately when idle instead
        # of spinning. Eliminates CPU burn in sequential sections between parallel
        # regions. May add ~1us wake-up latency per parallel region but avoids
        # hundreds of seconds of wasted spin time seen in profiling.
        VS_DEBUGGER_ENVIRONMENT "KMP_BLOCKTIME=0"
    )
endif()
