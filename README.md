# OpenPIC

OpenPIC is a 3D hybrid particle-in-cell plasma simulation framework with kinetic ions, fluid electrons, Lua-based problem configuration, and MPI/OpenMP parallel execution.

Ions are represented as kinetic macroparticles. Electrons are modeled as a massless fluid governed by a generalized Ohm's law. Electromagnetic fields are advanced on a staggered Yee C-grid using CGS-Gaussian units.

## Project Status

OpenPIC has been tested on Windows 11 and Ubuntu 26.04 LTS.

OpenPIC is research-oriented scientific software intended for plasma simulation, numerical experimentation, verification studies, and reproducible computational research.

The Windows build using MSVC 2022 is the primary verified configuration, with MPI support through Microsoft MPI. The Ubuntu build (GCC, OpenMPI) also works — see [Quick Start (Ubuntu)](#quick-start-ubuntu) — but has not been verified as thoroughly as Windows.

## Quick Start

The commands below are for Windows and should be run from a **Developer Command Prompt for Visual Studio 2022**. For Ubuntu/Linux, see [Quick Start (Ubuntu)](#quick-start-ubuntu) below.

### 1. Clone the Repository

```text
git clone https://github.com/dozzes/open-pic-3.git
cd open-pic-3
```

### 2. Build in Release Mode

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

The executable is generated at:

```text
build\bin\Release\open-pic.exe
```

### 3. Run a Verification Case

```text
build\bin\Release\open-pic.exe sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Diagnostic files are written to the `diag/` directory inside the simulation run directory.

For MPI builds, Linux notes, configuration options, and additional run modes, see the sections below.

## Quick Start (Ubuntu)

### 1. Install Build Tools

```bash
sudo apt update
sudo apt install -y build-essential cmake
```

### 2. Clone the Repository

```bash
git clone https://github.com/dozzes/open-pic-3.git
cd open-pic-3
```

### 3. Build in Release Mode

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

The executable is generated at:

```text
build/bin/open-pic
```

### 4. Run a Verification Case

Run with the path relative to the repository root:

```bash
./build/bin/open-pic sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

Or `cd` into the case directory and pass just the file name (the Lua path is
relative to the current directory, not the repo root):

```bash
cd sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM
~/open-pic-3/build/bin/open-pic main.lua
```

For an MPI-enabled build, install OpenMPI and reconfigure with
`-DOPENPIC_ENABLE_MPI=ON`:

```bash
sudo apt install -y libopenmpi-dev openmpi-bin
cmake -S . -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release -j
```

See [Linux Build](docs/OpenPIC_User_Guide.md#linux-ubuntu-build) in the User Guide for details.

## Key Features

* Three-dimensional hybrid particle-in-cell simulation
* Kinetic ion macroparticles
* Massless fluid-electron model
* Generalized Ohm's law
* Staggered Yee C-grid
* FDTD and PSTD field solvers
* Lua-based simulation configuration
* MPI distributed-memory parallelism
* OpenMP shared-memory parallelism
* Diagnostic output for particles and grid quantities
* Verification and regression-test workflows
* Separation of the numerical core from problem-specific simulation setup

## Physics Model

OpenPIC separates the computational core from problem definitions.

The numerical core is implemented in C++ and provides:

* particle initialization and particle pushing;
* charge and current deposition;
* electromagnetic-field solvers;
* particle-to-grid and grid-to-particle interpolation;
* boundary-condition infrastructure;
* diagnostic output;
* Lua bindings;
* MPI and OpenMP execution support.

Simulation scenarios are defined through Lua scripts. Geometry, initial conditions, boundary conditions, plasma parameters, particle populations, numerical methods, and diagnostic settings can therefore be changed without modifying the C++ core.

See [docs/OpenPIC_User_Guide.md](docs/OpenPIC_User_Guide.md) for the complete configuration reference.

One example problem supported by this architecture is the expansion of a plasma cloud into a magnetized background plasma in the transition between magnetolaminar and gas-dynamic regimes.

This problem can be characterized by two dimensionless parameters:

* **M_A = V_cloud / V_A** — Alfvén Mach number
* **Δ = (R_g / R_L)²** — magnetolaminar parameter

where:

* `V_cloud` is the characteristic cloud expansion velocity;
* `V_A` is the Alfvén velocity of the background plasma;
* `R_g` is the characteristic cloud radius;
* `R_L` is the characteristic ion Larmor radius.

The limiting regimes are:

* **Δ >> 1** — electromagnetic braking dominates;
* **Δ << 1** — gas-dynamic behavior dominates.

## Clone

```text
git clone https://github.com/dozzes/open-pic-3.git
cd open-pic-3
```

## Build

### Requirements

* CMake 3.13 or later
* C++17-compatible compiler

Windows:

* Microsoft Visual Studio 2022 with MSVC
* Microsoft MPI for MPI-enabled builds
* Run all commands from a **Developer Command Prompt for Visual Studio 2022**

Ubuntu/Linux:

* `build-essential` and `cmake` (`sudo apt install -y build-essential cmake`)
* `libopenmpi-dev` and `openmpi-bin` for MPI-enabled builds
* Run all commands from a regular terminal; see [Quick Start (Ubuntu)](#quick-start-ubuntu) for the full sequence

### Standard Build

Windows:

```text
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

The executable is generated at:

```text
build\bin\Release\open-pic.exe
```

Ubuntu/Linux:

```bash
sudo apt update
sudo apt install -y build-essential cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

The executable is generated at:

```text
build/bin/open-pic
```

### Build with MPI Support

Windows — install Microsoft MPI first:

```text
powershell -ExecutionPolicy Bypass -File tools\setup_msmpi.ps1
```

Ubuntu/Linux — install OpenMPI first:

```bash
sudo apt install -y libopenmpi-dev openmpi-bin
```

Configure and build OpenPIC with MPI enabled (same command on both platforms):

```text
cmake -S . -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release -j
```

See [docs/OpenPIC_User_Guide.md](docs/OpenPIC_User_Guide.md) for complete build instructions and troubleshooting information.

## Run

Run the uniform-equilibrium verification case.

Windows:

```text
build\bin\Release\open-pic.exe sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Ubuntu/Linux:

```bash
./build/bin/open-pic sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

The Lua path is relative to the current directory, not the repo root — if
you `cd` into the case directory first, pass just `main.lua`.

### MPI Execution

The following example starts four MPI ranks with four OpenMP threads per rank.

Windows:

```text
set OMP_NUM_THREADS=4
mpiexec -n 4 build\bin\Release\open-pic.exe -mpi sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Ubuntu/Linux:

```bash
OMP_NUM_THREADS=4 mpiexec -n 4 build/bin/open-pic -mpi sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

A successful run creates diagnostic output under the simulation run directory.

Diagnostics are written to the `diag/` directory as tab-separated `.dat` text files. Depending on the selected diagnostics, each row represents a grid node, particle, time sample, or other simulation quantity.

See [Diagnostics](docs/OpenPIC_User_Guide.md#11-diagnostics) in the User Guide for file formats and column definitions.

## Configuration

Simulations are configured through a Lua entry file, normally named `main.lua`.

The Lua configuration controls:

* grid dimensions and spacing;
* simulation time step;
* number of time steps;
* background plasma density;
* background electromagnetic fields;
* particle-group initialization;
* particle and field boundary conditions;
* FDTD or PSTD solver selection;
* diagnostic frequency and output selection.

Minimal example from:

`sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua`

```lua
local verify = require("verification_common")

verify.run({
    case_name = "UNIFORM_EQUILIBRIUM",
    nodes = 17,
    h = 1.0,
    Bz0 = 1.0,
    backgr_dens = 1.0,
    time_steps = 4,
    save_time_steps = 1,
    save_whole_grid = true,
})
```

See [Lua Configuration](docs/OpenPIC_User_Guide.md#6-lua-configuration) in the User Guide for the complete parameter reference.

## Verification and Testing

OpenPIC includes C++ tests, simulation-level verification cases, and regression-testing workflows.

The verification suite is intended to detect errors in:

* particle motion;
* field advancement;
* particle-grid coupling;
* equilibrium preservation;
* boundary-condition handling;
* MPI execution;
* diagnostic output;
* changes in numerical results between revisions.

Available verification and benchmark scenarios are located under:

```text
sim/tasks/04_VERIFICATION/
```

To build and run the C++ test suite (same commands on Windows and Linux):

```text
cmake -S . -B build-tests -DOPENPIC_BUILD_TESTS=ON
cmake --build build-tests --config Release --target openpic_tests
ctest --test-dir build-tests -C Release --output-on-failure
```

See:

* [Build and run C++ verification tests](docs/OpenPIC_User_Guide.md#build-and-run-c-verification-tests)
* [Regression testing](docs/OpenPIC_User_Guide.md#regression-testing)

for the complete verification and regression-checking workflow.

| Case | Directory | What it checks |
| --- | --- | --- |
| Uniform equilibrium (CIC) | `UNIFORM_EQUILIBRIUM` | A spatially uniform background plasma stays uniform: no spurious force from particle push, deposition, or field solve (default CIC particle shape) |
| Uniform equilibrium (NGP) | `UNIFORM_EQUILIBRIUM_NGP` | Same equilibrium check with nearest-grid-point deposition |
| Uniform equilibrium (TSC) | `UNIFORM_EQUILIBRIUM_TSC` | Same equilibrium check with triangular-shaped-cloud deposition |
| Linear Alfvén wave | `LINEAR_ALFVEN` | Numerical Alfvén-wave dispersion vs. the analytic Alfvén speed |
| Linear magnetosonic wave | `LINEAR_MAGNETOSONIC` | Numerical fast-magnetosonic dispersion vs. the analytic prediction |
| Field-loop advection | `FIELD_LOOP_DIVB` | div(B) = 0 preservation under the FDTD curl update while a magnetic field loop advects through the box |
| Orszag-Tang vortex | `ORSZAG_TANG_2D` | Qualitative reproduction of the classic 2D MHD turbulence benchmark |
| Whistler dispersion (coarse) | `WHISTLER_DISPERSION/COARSE` | Parallel whistler-branch ω(k) vs. the analytic formula at h = 2·d_i |
| Whistler dispersion (fine) | `WHISTLER_DISPERSION/FINE` | Same comparison at h = d_i, to check resolution convergence against COARSE |
| Hybrid resolution condition (valid) | `HYBRID_CONDITION_H_GT_CWPI/VALID` | Behavior when the hybrid-code condition h > c/ω_pi is satisfied |
| Hybrid resolution condition (violated) | `HYBRID_CONDITION_H_GT_CWPI/VIOLATED` | Contrast case with h < c/ω_pi |
| Inflow density conservation | `INFLOW_DENSITY_CONSERVATION` | Background density stays constant under continuous inflow/outflow particle injection |
| Resonant ion-beam instability | `ION_BEAM_RESONANT` | Growth rate of the electromagnetic ion/ion resonant instability vs. its closed-form analytic estimate |

Each case is run and compared against its analytic reference manually (see the case's `main.lua` header for the derivation); there is no automated pass/fail gate or stored benchmark output in the repository yet.

## Example Results

Add one or more representative OpenPIC results here, such as:

* plasma-density slices;
* magnetic-field distributions;
* expanding-cloud fronts;
* particle distributions;
* conservation-error plots;
* comparisons with analytical or published results;
* MPI scaling measurements.

Example Markdown:

```markdown
![Plasma density distribution](docs/images/plasma-density.png)

*Example plasma-density distribution produced by OpenPIC.*
```

## Repository Layout

```text
src/                  C++ source files and headers
tests/                Lightweight C++ verification test executable (openpic_tests)
sim/
  lib/                Shared Lua simulation library
                      (simulation core, callbacks, boundary conditions)
  tasks/
    04_VERIFICATION/  Verification and benchmark cases
docs/
  OpenPIC_User_Guide.md
                      Complete user guide
tools/                Setup scripts and diagnostic utilities
3rdparty/             Vendored third-party libraries
                      (Lua, sol2, fmt, pocketfft)
```

## Third-Party Dependencies

| Library                                            | Version | Purpose                                    |
| -------------------------------------------------- | ------: | ------------------------------------------ |
| [Lua](https://www.lua.org/)                        |   5.4.8 | Simulation configuration and scripting     |
| [sol2](https://github.com/ThePhD/sol2)             |   3.3.0 | C++ and Lua bindings                       |
| [fmt](https://github.com/fmtlib/fmt)               |  10.2.1 | String formatting                          |
| [pocketfft](https://github.com/mreineck/pocketfft) |       — | FFT implementation used by the PSTD solver |

## Documentation

* [OpenPIC User Guide](docs/OpenPIC_User_Guide.md) — build instructions, Lua configuration, boundary conditions, MPI execution, diagnostics, testing, and troubleshooting
* [Shared Lua Library Reference](sim/lib/README.md) — reusable Lua simulation components and APIs

## Citation

> **3D modeling of the plasma expansion in a dipole magnetic field.**  
> *Proceedings of the XXX International Conference on Phenomena in Ionized Gases (ICPIG 2011)*, Belfast, Northern Ireland, UK, August 2011.

## Authors and Acknowledgments

OpenPIC was developed at the Institute of Radiophysics and Electronics of the National Academy of Sciences of Armenia.

Lead developer:

* David Osipyan

## License

OpenPIC is distributed under the terms of the [BSD 3-Clause License](LICENSE).
