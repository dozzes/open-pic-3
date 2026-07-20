# OpenPIC

3D hybrid particle-in-cell code for modeling plasma clouds in a magnetized background.

Ions are treated as macroparticles (kinetic). Electrons are a massless fluid governed by a generalized Ohm's law. Fields are advanced on a staggered Yee C-grid in CGS-Gaussian units.

## Physics

The code targets the magnetolaminar / gas-dynamic transition regime. The key dimensionless parameters are:

- **Ma** = V_cloud / V_A — Alfven Mach number
- **Delta** = (R_g / R_L)² — magnetolaminar parameter  
  Delta >> 1: electromagnetic braking dominates  
  Delta << 1: gas-dynamic regime

Available magnetic field solvers:
- **FDTD** — standard Yee finite-difference time-domain (default, production)
- **PSTD** — pseudo-spectral time-domain (numerical experiment, lower dispersion)

## Clone

This repository uses two submodules: `sim/tasks/examples` (example cases) and
`docs/research` (research notes, analysis, article drafts). Verification
cases and the user guide live directly in this repository.

Fresh clone:

```powershell
git clone --recurse-submodules https://github.com/dozzes/open-pic-2.git
```

If you already cloned without `--recurse-submodules`:

```powershell
git submodule update --init --recursive
```

Pulling later, in an existing clone:

```powershell
git pull
git submodule update --init --recursive
```

## Build

Requirements: CMake 3.20+, MSVC 2022 (Windows) or GCC/Clang (Linux), C++17.

```powershell
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Output: `build/bin/Release/open-pic.exe`

**With MPI support:**

```powershell
# Install MS-MPI first (Windows)
powershell -ExecutionPolicy Bypass -File tools\setup_msmpi.ps1

cmake -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release
```

See [docs/OpenPIC_User_Guide.md](docs/OpenPIC_User_Guide.md) for full build and run instructions.

## Run

```powershell
build\bin\Release\open-pic.exe sim\tasks\examples\01_BASELINE\MA2_D44_UNIF_COLD\main.lua
```

With MPI (4 ranks, 4 OpenMP threads each):

```powershell
$env:OMP_NUM_THREADS = 4
mpiexec -n 4 build\bin\Release\open-pic.exe -mpi sim\tasks\examples\01_BASELINE\MA2_D44_UNIF_COLD\main.lua
```

Diagnostics are written to `diag/` inside the run directory.

## Repository Layout

```
src/              C++ source and headers
sim/
  lib/            Shared Lua library (sim_core, callbacks, boundary conditions)
  tasks/
    04_VERIFICATION/  Verification and benchmark cases
    examples/         [submodule: openpic-examples] Example simulation cases
      01_BASELINE/    Reference cases: uniform field, Ma=2, Delta=44
      02_BG_FLOW/     Cases with flowing background plasma
      03_DIPOLE/      Dipole magnetic field geometry
      SOLAR_WIND/     Background-only wind reference (no cloud)
docs/
  OpenPIC_User_Guide.md  User guide
  research/              [submodule: openpic-docs] Research notes, analysis, article drafts
tools/            Setup scripts and diagnostic utilities
3rdparty/         Vendored libraries (Lua, sol2, fmt, pocketfft)
```

### Task naming convention

```
MA{N}_D{delta}_{FIELD}_{COLD|THERM}[_BG{M}][_modifiers]
```

| Segment | Meaning |
|---------|---------|
| `MA{N}` | Alfven Mach number |
| `D{delta}` | Delta parameter (P = decimal point, e.g. D0P8) |
| `UNIF` / `DIP` / `PLANET_DIP` | Field geometry |
| `COLD` / `THERM` | Electron model |
| `BG{M}` | Background flow Mach number (if non-zero) |

## Third-party dependencies

| Library | Version | Use |
|---------|---------|-----|
| [Lua](https://www.lua.org/) | 5.4.8 | Configuration scripting |
| [sol2](https://github.com/ThePhD/sol2) | 3.3.0 | C++/Lua binding |
| [fmt](https://github.com/fmtlib/fmt) | 10.2.1 | String formatting |
| [pocketfft](https://github.com/mreineck/pocketfft) | — | FFT for PSTD solver |

## Documentation

- [docs/OpenPIC_User_Guide.md](docs/OpenPIC_User_Guide.md) — full reference: build, Lua config, boundary conditions, MPI mode, diagnostics, troubleshooting
- [sim/tasks/README.md](sim/tasks/README.md) — task group descriptions and naming convention
- [sim/lib/README.md](sim/lib/README.md) — shared Lua library reference
- `docs/research/` — research notes, analysis, article drafts (submodule; see Clone section above)
