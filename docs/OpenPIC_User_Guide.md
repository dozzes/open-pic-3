# OpenPIC User Guide

OpenPIC has been tested on Windows 11 and Ubuntu 26.04 LTS.

## 1. Overview

OpenPIC is a hybrid particle-in-cell code for modeling ion kinetic dynamics in a magnetized plasma. Ions are represented by macroparticles, while electrons are treated as a massless fluid. The electromagnetic fields are advanced on a staggered C-grid, and the electric field is computed from a generalized Ohm's law.

The code is driven by Lua input files. A typical simulation directory contains the main physical setup, boundary-condition scripts, diagnostic settings, and a run script. The executable reads the Lua configuration, initializes particle groups and fields, advances the simulation, and writes grid and particle diagnostics.

## 2. Getting the Code

To checkout the OpenPIC source, do the following:

```text
git clone https://github.com/dozzes/open-pic-3.git
cd open-pic-3
```

## 3. Build

OpenPIC can be built and run on Windows 11 and Ubuntu 26.04.

OpenPIC uses the CMake build system. To configure a build, do the following from the top-level source directory:

Windows:

```text
cmake -S . -B build
cmake --build build --config Release
```

The Release executable is produced at:

```text
build\bin\Release\open-pic.exe
```

### Linux (Ubuntu) Build

Install a compiler toolchain and CMake:

```bash
sudo apt update
sudo apt install -y build-essential cmake
```

Configure and build from the top-level source directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

The Release executable is produced at:

```text
build/bin/open-pic
```

The Lua file path is relative to the current directory, not the repo root:

```bash
./build/bin/open-pic sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

Or run it from inside the case directory:

```bash
cd sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM
/path/to/open-pic-3/build/bin/open-pic main.lua
```

For an MPI-enabled Linux build, install OpenMPI first (or run
`tools/setup_mpi_linux.sh`), then reconfigure with `-DOPENPIC_ENABLE_MPI=ON`:

```bash
sudo apt install -y libopenmpi-dev openmpi-bin
cmake -S . -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release -j
```

See [§14, MPI Distributed Mode](#14-mpi-distributed-mode) for running MPI
builds and the tradeoffs of the replicated-grid model.

### CMake build options

All optional OpenPIC CMake switches are disabled by default, so the plain
Release build shown above (no extra `-D` flags) is the production build:

| Option | Default | Use for | Production impact |
|---|---:|---|---|
| `OPENPIC_BUILD_TESTS` | `OFF` | Builds `openpic_tests` and enables CTest verification tests. | No effect on `open-pic.exe`; adds a separate test executable. |
| `OPENPIC_ENABLE_MPI` | `OFF` | Builds MPI-capable OpenPIC. Required for `open-pic.exe -mpi ...`. | Adds MPI dependency and MPI code paths. Without `-mpi`, the executable still uses the single-process path. |
| `OPENPIC_PROFILE_PARTICLE_PHASES` | `OFF` | Prints detailed timing blocks for particle phases and high-level simulation phases. | Profiling only. Do not use for production timing except controlled performance tests. |
| `OPENPIC_PARTICLE_STATIC_SCHEDULE` | `OFF` | Replaces the particle half-step OpenMP loop schedule with `schedule(static)` for scaling experiments. | Can change performance characteristics. Treat as experimental until validated for production workloads. |
| `OPENPIC_NATIVE_ARCH` | `ON` | GCC/Clang builds only (no effect on MSVC, which always builds with `/arch:AVX2`): adds `-march=native`. | Maximizes throughput but makes the binary's numerical results machine-dependent. Disable for verification/validation runs that must reproduce bit-for-bit across hosts. |

MPI build (same command on Windows and Linux):

```text
cmake -S . -B build-mpi -DOPENPIC_ENABLE_MPI=ON
cmake --build build-mpi --config Release
```

OpenMP profiling build. Windows (`^` continues the line in cmd.exe):

```text
cmake -S . -B build-omp-profile ^
  -DOPENPIC_PROFILE_PARTICLE_PHASES=ON ^
  -DOPENPIC_PARTICLE_STATIC_SCHEDULE=ON ^
  -DOPENPIC_ENABLE_MPI=OFF
cmake --build build-omp-profile --config Release
```

Linux (`\` continues the line in the shell):

```bash
cmake -S . -B build-omp-profile \
  -DOPENPIC_PROFILE_PARTICLE_PHASES=ON \
  -DOPENPIC_PARTICLE_STATIC_SCHEDULE=ON \
  -DOPENPIC_ENABLE_MPI=OFF
cmake --build build-omp-profile --config Release
```

Do not mix profiling builds with production baselines. The profiling code is
compiled out unless `OPENPIC_PROFILE_PARTICLE_PHASES=ON`, but when enabled it
adds timers and log output. `OPENPIC_PARTICLE_STATIC_SCHEDULE=ON` is also a
deliberate performance experiment, not the default production path.

### Runtime OpenMP setting

`OMP_NUM_THREADS` is not a CMake option. It is a runtime environment variable
that controls how many OpenMP worker threads the executable uses.

Windows:

```text
set OMP_NUM_THREADS=16
build\bin\Release\open-pic.exe main.lua
```

Linux:

```bash
export OMP_NUM_THREADS=16
build/bin/open-pic main.lua
```

On a 24-thread workstation, `16` is the recommended setting. `24` is a bit
faster but uses the CPU less efficiently.

### Build and run C++ verification tests

OpenPIC has a lightweight C++ verification test executable that does not depend
on GoogleTest or other external test frameworks. It is disabled by default and
enable it with (same commands on Windows and Linux):

```text
cmake -S . -B build-tests -DOPENPIC_BUILD_TESTS=ON
cmake --build build-tests --config Release --target openpic_tests
ctest --test-dir build-tests -C Release --output-on-failure
```

The first test target is `openpic_tests`. It currently checks:

- **Particle push:** Boris pusher in a uniform electric field against the
  analytic velocity update; Boris pusher in a uniform magnetic field for
  speed conservation.
- **Deposition/gather:** CIC scatter weight conservation for `NP`, `UP`, and
  `UP_NP`; NGP and quadratic TSC gather/deposit consistency on the staggered
  grid; partition-of-unity of the gather weights; face-velocity scatter.
- **Field solvers:** FDTD curl dispersion symbol; the Faraday update
  preserves div(B) = 0; PSTD plane-wave Faraday derivative against the
  analytic spectral derivative.
- **Linear algebra / FFT:** Thomas tridiagonal solve vs. dense solve;
  pocketfft complex-to-complex round trip.
- **Infrastructure:** grid index/centering conventions (`to_idx`); RNG range
  utility; MPI serial-mode defaults and no-op density reduction; a full
  single-process `simulate()` smoke run.

### Run verification benchmark tasks

Each case is a self-contained run directory, run individually with the built
executable.

Windows:

```text
build\bin\Release\open-pic.exe sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Linux:

```bash
build/bin/open-pic sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

Then inspect the resulting `diag/*.dat` files in that task's own directory
(see [Diagnostics](#11-diagnostics)) against the expectation described for
that case below. All 11 cases are single-process only (`verify.run()`
enforces this) and most complete in seconds; `ION_BEAM_RESONANT` is the
exception (tens of thousands of steps, on the order of hours).

### Verification tasks

What each case under `sim/tasks/04_VERIFICATION/` checks, from its own
`main.lua`:

- **`UNIFORM_EQUILIBRIUM`** (and the **`UNIFORM_EQUILIBRIUM_NGP`** /
  **`UNIFORM_EQUILIBRIUM_TSC`** variants): a spatially uniform background
  plasma in a uniform `Bz0` should stay uniform after a few push/deposit/
  field-solve cycles (4 steps) — a basic equilibrium/no-self-excitation
  check. The three variants exercise the three particle-shape functions:
  standard CIC, `ScatterAlg.NGP`, and `ScatterAlg.TSC`.
- **`FIELD_LOOP_DIVB`**: seeds a localized circular in-plane `B`
  perturbation (a "field loop") on top of a uniform background with a
  small flow, over 8 steps — checks that the FDTD Faraday update does not
  introduce spurious field structure from a compact, non-trivial `B` seed.
- **`LINEAR_ALFVEN`**: seeds a small-amplitude sinusoidal transverse
  velocity/`B` perturbation aligned with `B0` (a linear parallel-Alfven-wave
  initial condition), 8 steps — sanity check of linear, non-compressive
  wave propagation.
- **`LINEAR_MAGNETOSONIC`**: seeds a compressive sinusoidal density/`|B|`
  perturbation transverse to the propagation direction (a linear
  magnetosonic-wave-like initial condition), 8 steps — sanity check of
  linear, compressive wave propagation.
- **`ORSZAG_TANG_2D`**: the classic Orszag-Tang vortex initial condition
  (sinusoidal velocity and magnetic field in the x-y plane), run for 8
  steps — a short smoke check that this well-known non-trivial 2D
  initial condition evolves without immediately blowing up; not a full
  nonlinear-turbulence validation at this step count.
- **`HYBRID_CONDITION_H_GT_CWPI/VALID`** and **`.../VIOLATED`**: the same
  parallel-Alfven-wave seed as `LINEAR_ALFVEN`, run at two background
  densities so that the grid spacing `h` is either larger than the ion
  inertial length `c/omega_pi` (`VALID`, printed as `h/(c/Wpi)`) or smaller
  than it (`VIOLATED`). Documents/exercises the hybrid-PIC resolution
  validity condition `h > c/omega_pi` rather than asserting a pass/fail
  physics result by itself.
- **`INFLOW_DENSITY_CONSERVATION`**: audits the X-axis particle-injection
  mechanism (`inject_background_particles` / `add_inflow_plane_x`). A
  uniform background stream drifts through a field-free box (`Bz0 = 0`)
  and exits at the far wall; the check is whether inflow replenishment
  keeps the interior density at `backgr_dens` over 2000 steps (5 full
  box-crossings) instead of draining. Long relative to the other cases
  because it needs several steady-state injection cycles, not one wave
  period.
- **`ION_BEAM_RESONANT`**: a two-population setup — a drifting ion "beam"
  as the `backgr` group (using the x-inflow mechanism above to stay
  replenished) plus a static "core" population as an `extra_group` — with
  a circularly-polarized seed field tuned to the resonance condition of
  the electromagnetic ion/ion beam instability. This is the hybrid-PIC
  analog of the two-stream instability (OpenPIC's electrons are an exactly
  massless fluid, so no electrostatic two-stream mode exists here). Checks
  for genuine exponential growth of the transverse field at the predicted
  resonant wavenumber with the correct polarization sense, over ~47000
  steps (~5 e-folds) — by far the most expensive of the 11 cases.
- **`WHISTLER_DISPERSION/COARSE`** and **`.../FINE`**: seeds a
  circularly-polarized `B` perturbation around `Bz0` (selecting only the
  right-hand/whistler branch) and compares the numerically measured
  rotation period (from `Bx` zero-crossings in the saved grid time series)
  against the analytic parallel whistler dispersion relation `omega = V_A
  k (sqrt(1+(k d_i/2)^2) + k d_i/2)`. `COARSE` uses `h = 2 d_i`, `FINE`
  uses `h = d_i` (same `k`, same box) — comparing the two also checks
  whether finer resolution improves dispersion accuracy in the regime
  where the hybrid approximation is least valid.

None of these cases currently print an automated `PASS`/`FAIL` verdict;
checking the result means reading the relevant `diag/*_grd_*.dat` time
series (as described per case above) against the expectation stated in
that case's own `main.lua` header comment. For diffing a run's output
against a stored baseline instead of an analytic expectation, see
[Regression testing](#regression-testing) below.

### Regression testing

Besides `ctest` (unit-level, checks physics kernels in isolation) and the
verification benchmark suite above (checks whole runs against analytic or
qualitative expectations), OpenPIC ships tools that diff the actual grid
output of two runs against each other. This is the tool to use when you want
to confirm that a code change did not silently alter results — for example,
after a refactor that is supposed to be behavior-preserving.

Compare all `*grd*.dat` files under two run directories, cell by cell, within
an absolute/relative tolerance.

Windows:

```text
powershell -ExecutionPolicy Bypass -File tools\Compare-GridOutputs.ps1 -RunA path\to\run_a -RunB path\to\run_b
```

`tools/compare_grid_outputs.py` is a Python port of the same tool, for Linux,
macOS, or Windows with Python installed:

```bash
python tools/compare_grid_outputs.py path/to/run_a path/to/run_b --atol 1e-10 --rtol 1e-8
```

Both report per-file `PASS`/`FAIL` with the maximum absolute/relative
deviation and its location, and exit non-zero on any mismatch, so they can be
used as a CI gate.

`tools\Compare-BaselineAfterFix.ps1` wraps the same comparison around a
before/after workflow: it backs up the existing `diag/` output of a case,
re-runs the case with the freshly built executable, diffs the new output
against the backup, and reports the result. It defaults to a case path that
lives in the (non-public) examples set — pass `-CaseDir` to point it at a
case that exists in this repository, e.g.
`sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM`.

For MPI-specific regression checks (single-process output vs. MPI output),
see `tools\Run-MpiTaskGridChecks.ps1` in [§14, Regression Testing
MPI](#regression-testing-mpi).

## 4. Hardware Requirements

OpenPIC keeps the full field grid and the full macroparticle set in memory
(no out-of-core or streaming I/O), so RAM is the resource most likely to
limit problem size. Two structures dominate:

- **`Particle`** (one macroparticle): position, velocity, weight, group id —
  cache-line aligned to 64 bytes. Memory for particles scales as
  `N_particles x 64 bytes`, independent of grid size or thread count.
- **grid node** (`Cell`): density, `B`, `E`, `UE`, `UP`, and their face/edge
  companions — roughly 150 bytes per node (exact value depends on compiler
  padding). Memory for the grid scales as `N_cells x ~150 bytes`.

OpenMP scatter uses one thread-private density accumulator (`~56 bytes` per
node) per worker thread to avoid write races, so grid-related memory also
scales with `OMP_NUM_THREADS`:

```text
grid_bytes ~= N_cells * (150 + 56 * OMP_NUM_THREADS)
particle_bytes ~= N_particles * 64
total_bytes ~= grid_bytes + particle_bytes
```

Worked example: a 100x100x100 grid (1e6 cells) with `OMP_NUM_THREADS=16` and
10 million macroparticles:

```text
grid:      1e6 * (150 + 56*16)  ~= 1.05 GB
particles: 1e7 * 64             ~= 0.64 GB
total                           ~= 1.7 GB
```

This is a lower bound on the working set (core simulation state only; it
excludes the executable, Lua runtime, OS overhead, and diagnostic I/O
buffers) — treat it as an order-of-magnitude planning figure, not an exact
budget. In [MPI mode](#14-mpi-distributed-mode) the grid is replicated on
every rank (multiply `grid_bytes` by rank count), while particle memory is
divided across ranks.

For thread-count-vs-wall-time trade-offs (as opposed to memory), see
[OpenMP Performance Testing](#13-openmp-performance-testing).

## 5. Run Directory Layout

A run directory usually contains:

```text
main.lua          Main physical and numerical setup
bound_cond.lua    Boundary-condition callbacks
em_field.lua      Initial electromagnetic field setup
parts_count.lua   Particle count helpers
print_mpi.lua     MPI-aware printing helpers
utils.lua         Utility functions
run.bat           Local launch script
diag/             Grid diagnostics
```

The code writes diagnostic files such as:

```text
diag/all_z_27_grd_0000.dat
diag/backgr_z_27_grd_0000.dat
diag/cloud_z_27_grd_0000.dat
```

`all_*` diagnostics contain moments of all particles. Group-specific files, such as `backgr_*` and `cloud_*`, contain moments for individual particle groups.

## 6. Lua Configuration

OpenPIC simulations are configured entirely in Lua; there is no other input
format. Two workflows exist:

- **Production workflow** (`require("sim_core").run()`): used for physical
  scenario cases such as the cloud-in-background-plasma template in
  `sim/lib/main.lua`. Grid spacing `h` and time step `tau` are *derived*
  from physical resolution parameters (Larmor radius, CFL/cyclotron
  stability limits) rather than set directly.
- **Verification workflow** (`require("verification_common").run()`, in
  `sim/tasks/04_VERIFICATION/lib/verification_common.lua`): used for the
  test cases under `sim/tasks/04_VERIFICATION/`. Grid spacing `h` is a
  direct input parameter; `tau` is auto-computed from stability limits
  unless overridden via `params.tau`. This workflow is single-process only.

The groups below use the parameter names from `sim/lib/main.lua`
(production workflow); the verification workflow's `verify.run({...})`
table uses the same names for the parameters it exposes — see the minimal
example in [Configuration](../README.md#configuration).

### Case label

```lua
case_name = "Case A: cold electrons, stationary background, FDTD solver"
```

Human-readable label used to identify output; has no effect on the physics.

### Background plasma

```lua
backgr_mi   = 1 * cgs.mp   -- background ion mass
backgr_dens = 1.0e+15      -- background ion number density [cm^-3]
backgr_Z    = 1            -- background ion charge number

backgr_flow_Ma_x = 0.0     -- background drift velocity, in units of V_A0
backgr_flow_Ma_y = 0.0
backgr_flow_Ma_z = 0.0
```

All-zero flow gives a stationary background. A nonzero `backgr_flow_Ma_x`
gives a flowing background (e.g. `8` for a solar-wind-like super-Alfvenic
flow) — see [Stationary and Flowing
Backgrounds](#10-stationary-and-flowing-backgrounds). Particle injection at
the inflow boundary is currently implemented only for X-directed flow.

### External magnetic field

```lua
Bz0 = 100.0   -- background magnetic field [Gauss], along Z
```

The Alfven speed `V_A0` used throughout the case (Mach number
normalization, resolution formulas) is derived from `Bz0` and
`backgr_dens`.

### Plasma cloud

```lua
cloud_mi       = 1 * cgs.mp   -- cloud ion mass
cloud_ions_num = 6.0e+18      -- total number of real (physical) ions in the cloud
cloud_Z        = 1            -- cloud ion charge number
Ma             = 2.0          -- Alfven Mach number: cloud_v_max / V_A0
```

This block is specific to the cloud-expansion example scenario; other Lua
scenarios can omit it and populate particle groups differently — see
[Particle Groups](#8-particle-groups).

### Algorithm selection

```lua
magnetic_field_alg = MagneticFieldAlg.FDTD   -- or MagneticFieldAlg.PSTD
cold_electrons_enabled = true                -- false enables Spitzer resistivity + heat conduction
cloud_jitter_enabled  = false                -- randomize macroparticle lattice positions
backgr_jitter_enabled = false
pstd_quick_test_steps = nil                  -- cap step count for a PSTD sanity check; nil = full run
```

See [PSTD status](#pstd-status) for FDTD-vs-PSTD guidance, [Electron
Model](#7-electron-model) for the cold/finite-temperature electron closure,
and [Console verbosity](#console-verbosity) for
`pic_config.set_verbosity_level()`.

### Resolution parameters

```lua
points_on_RL = 10   -- grid cells per ion Larmor radius (minimum ~5, production ~10)
h_R0 = 4.0           -- initial cloud radius, in grid steps
break_times = 4.0    -- simulated duration, in magnetic braking times T_b = R_b / V_max

cloud_parts_on_step  = 25   -- macroparticles per grid step per axis, cloud
backgr_parts_on_step = 4    -- macroparticles per grid step per axis, background
```

In this workflow, `h` and `tau` are not set directly — they are derived
from these resolution parameters plus the CFL and cyclotron stability
limits (see `sim/lib/physics.lua`). The verification workflow instead
takes `h` as a direct parameter (`nodes`/`h` in `verify.run({...})`) with
`tau` auto-computed unless overridden.

Particle count scales as `parts_on_step^3` times the populated volume, so
it grows quickly — this is the main memory/cost knob; see [Hardware
Requirements](#4-hardware-requirements).

### Launching the run

```lua
require("sim_core").run()
```

Calls, in order: physics derivation (resolution -> `h`/`tau`, normalization
constants) -> particle initialization -> field initialization -> the main
time-stepping loop -> diagnostic output.

### Electron closure overrides

Three electron-fluid closure coefficients — resistivity `eta`, heat
conduction `chi`, interspecies friction `kappa` — have physically-motivated
defaults assigned inside `run()`. Overrides only take effect if set AFTER
the `require("sim_core").run()` call:

```lua
pic_parameters.Spitzer_eta_floor_mult = 1.0             -- floor eta at base_eta above Te_ref
pic_parameters.Spitzer_Te_ref = 0.0                      -- diagnostic: constant eta everywhere
pic_parameters.electron_thermal_conductivity = 1.93e+5   -- physical Spitzer chi at Te=1eV, n=1e15 cm^-3
pic_parameters.kappa_friction = 0.1                       -- legacy interspecies-friction closure
pic_parameters.kappa_B_ref    = Bz0
```

`kappa_friction`/`kappa_B_ref` together replace the Spitzer `eta` model
entirely, rather than composing with it. Full derivations and measured
effects: `sim/lib/README.md`, section "Electron closure parameters".

### Numerical parameters

```lua
h                         -- grid spacing
tau                       -- time step
pic_parameters.time_steps
pic_parameters.save_time_steps
pic_parameters.dens_cutoff
```

Particle pushing and scattering are selected through Lua-bound enums:

```lua
pic_parameters.push_method = ParticlePushAlg.Boris
pic_parameters.scatter_method = ScatterAlg.Standard
```

`ScatterAlg.Standard` is the historical CIC particle shape. For controlled
hybrid cancellation-error studies, select NGP consistently for deposition and
field gathering through the case runner:

```lua
require("sim_core").run({ scatter_alg = ScatterAlg.NGP })
```

For a smoother quadratic particle shape with three-point support, use:

```lua
require("sim_core").run({ scatter_alg = ScatterAlg.TSC })
```

The selection applies separately on each staggered C-grid component and uses
the same shape for deposit and gather; CIC remains the default for all
existing cases. A controlled translation relative to the grid can be combined
with any shape:

```lua
require("sim_core").run({
    scatter_alg = ScatterAlg.TSC,
    cloud_center_shift_h = { x = 0.25, y = 0.25, z = 0.0 },
})
```

The magnetic-field update can be selected as:

```lua
magnetic_field_alg = MagneticFieldAlg.FDTD
-- or
magnetic_field_alg = MagneticFieldAlg.PSTD
```

Use FDTD as the reference production method unless a PSTD-specific comparison is being performed.

### PSTD status

`PSTD` is currently an experimental/research solver, not the default production
path. It is useful for controlled comparisons against FDTD numerical dispersion
because it replaces the finite-difference curl with a spectral derivative.

Important limitations:

- PSTD uses global FFTs and is therefore non-local.
- PSTD usually requires a smaller time step (`PSTD_tau`) than FDTD production
  tasks.
- In MPI replicated-grid mode, every rank performs the full FFT on the full
  grid, so PSTD is not a good production MPI solver.
- PSTD has kernel-level verification coverage, but not yet a full production
  benchmark suite.

Current automated guard:

```text
openpic_tests: pstd_plane_wave_Faraday_derivative
```

This unit test initializes a sinusoidal `Ez(y + h/2)` mode and checks that the
PSTD Faraday half-step updates `Bx` with the analytic spectral derivative.
There is no Lua-level PSTD baseline smoke case in this repository; PSTD
coverage is currently limited to the kernel-level unit test above, so do not
treat PSTD as validated at the whole-simulation level yet.

### Console verbosity

`pic_config.set_verbosity_level(level)`, called before `require("sim_core").run()`,
controls per-step console output only — it has no effect on diagnostic files.
Only MPI rank 0 prints.

| Level | Output |
|-------|--------|
| **0** | Silent |
| **1** (default) | A single-line progress bar (percent, elapsed, ETA, step/s), redrawn in place on a terminal |
| **2** | Level 1, plus a `print_tm()` timer line after each named phase (field solve, particle push, etc.) — for profiling |

Level 1:

```text
  [===========================>]  94.74%  270/285  elapsed 08:15  ETA 00:27  0.5 step/s
```

Level 2 additionally prints lines like this between progress updates:

```text
 366 sec: move_particles_half_time
```

## 7. Electron Model

The electron model can be configured as a cold-electron reference or as a finite-temperature electron fluid.

Cold-electron reference:

```lua
cold_electrons_enabled = true
pic_parameters.resistivity = 0.0
pic_parameters.electron_thermal_conductivity = 0.0
```

Finite-temperature electrons:

```lua
cold_electrons_enabled = false
pic_parameters.electron_thermal_conductivity = 1e-12 * (h*h) / tau
```

For thermal transport scans, useful values are:

```lua
0.0
1e-12 * (h*h) / tau
1e-11 * (h*h) / tau
5e-11 * (h*h) / tau
```

The heat-conductivity coefficient should be treated as a controlled model parameter and kept fixed within a single production case.

## 8. Particle Groups

Particle groups are created in Lua:

```lua
pic_particle_groups:create_group("backgr", 1.0, 1.0, Diagnostics.save_grid_values)
pic_particle_groups:create_group("cloud",  1.0, 1.0, Diagnostics.save_grid_values)
```

A moving or reservoir-like background can be marked explicitly:

```lua
pic_particle_groups:set_boundary_kind("backgr", BoundaryKind.FlowBackground)
```

This lets boundary-condition callbacks treat the background group differently from ordinary particle groups without hard-coding group names in C++.

## 9. Boundary Conditions

Boundary conditions are implemented as Lua callbacks called from C++:

```lua
on_set_boundary_NP()
on_set_boundary_UP()
on_set_boundary_UE()
on_set_boundary_EF()
on_set_boundary_MF()
on_set_boundary_Te()
on_set_boundary_eta()
```

Group-specific callbacks receive the group name:

```lua
on_set_boundary_group_NP(density_grid, group_name)
on_set_boundary_group_UP(density_grid, group_name)
```

Important convention:

- `NP`, `Te`, and `eta` are cell-centered.
- `Ex`, `UEx`, and `UPx` are X-face-centered.
- `Ey`, `UEy`, and `UPy` are Y-face-centered.
- `Ez`, `UEz`, and `UPz` are Z-face-centered.
- `Bx`, `By`, and `Bz` live on grid edges.

Lua boundary conditions should set physical values. In particular, `NP` should be set as a physical density, not as `density * cell_volume`. The C++ finalization pipeline normalizes particle weights before calling the final `NP` boundary condition.

## 10. Stationary and Flowing Backgrounds

A stationary background is selected by:

```lua
backgr_flow_Ma_x = 0.0
backgr_flow_Ma_y = 0.0
backgr_flow_Ma_z = 0.0
```

A flowing background along X is selected by:

```lua
backgr_flow_Ma_x = 0.1
backgr_flow_Ma_y = 0.0
backgr_flow_Ma_z = 0.0
```

For a flowing background, background particles must be replenished at the inflow boundary. In the current Lua setup, particle injection is implemented for X-directed flow. If all background velocity components are zero, the injection function returns immediately and no new particles are added.

## 11. Diagnostics

All diagnostic output is plain tab-separated text (`std::ofstream`/`fwrite`,
one header row plus one row per grid node or particle) — there is no
binary, HDF5, VTK, or comma-separated output. Files are written directly
into `diag/` inside the run directory, one file per group/plane/snapshot;
there is no single combined output file.

### Dump frequency

A dump is written on time step `n` whenever `n % save_time_steps == 0`, plus
always on the final step (`n == time_steps`). `save_time_steps` is a
top-level Lua parameter (`pic_parameters.save_time_steps` /
`save_time_steps` in `verify.run({...})`); see [Numerical
parameters](#numerical-parameters).

### Grid files

Naming: `{group}_grd_{step:0N}.dat` for a whole-grid dump
(`save_whole_grid = true`), or `{group}_{x|y|z}_{level}_grd_{step:0N}.dat`
for a single-plane dump (`save_grid_x_plains` / `_y_plains` / `_z_plains`),
e.g. `all_z_27_grd_0000.dat` is the Z=27 plane of the `all` group at step
0. `{step:0N}` is zero-padded to the width of `time_steps`. `{group}` is
`all` for the whole-simulation moments, or a particle group name (e.g.
`backgr`, `cloud`) for per-group moments — see [Particle
Groups](#8-particle-groups).

Column layout (one header row, then one row per grid node):

```text
X Y Z NP B Bx By Bz E Ex Ey Ez div(E) UP UPx UPy UPz Ma UE UEx UEy UEz Te c_Wpi
```

Typical quantities:

- `NP`: ion density normalized by the reference density.
- `B`, `Bx`, `By`, `Bz`: magnetic field normalized by the reference field.
- `E`, `Ex`, `Ey`, `Ez`: electric field normalized by the reference field scale.
- `div(E)`: diagnostic normalized as `div(E)/(E_scale/L_scale)`.
- `UP`: ion fluid velocity.
- `UE`: electron fluid velocity.
- `Te`: electron temperature normalized by the reference electron temperature.
- `c_Wpi`: ion inertial length estimate.

For boundary checks, inspect the first and last saved grid layers in `X` and `Y`. In a stationary-background cold-electron run, edge `NP` should remain close to unity and edge `UPx` should remain close to zero.

### Particle files

Naming: `{group}_parts_{step:0N}.dat`, written for any particle group whose
diagnostics flags include `save_positions`, plus `all_parts_{step:0N}.dat`
when `save_all_particles = true`. Unlike the grid files, particle files
have **no header row** — columns are, in order:

```text
X Y Z V Vx Vy Vz B Bx By Bz E Ex Ey Ez
```

(position, speed and velocity components, local field magnitude and
components at the particle's position). Inactive (removed/absorbed)
particles are skipped, so row count can be less than the particle count at
`t=0`.

Marker-particle files (`{marker_set}_markers_{step:0N}.dat`, for particles
tracked individually via `MarkerParticles`) use the same column set but
**do** include a header row.

### Tools

- `tools\Compare-GridOutputs.ps1` / `tools\compare_grid_outputs.py` — diff
  grid files between two runs; see [Regression testing](#regression-testing).
- `tools\Compare-BaselineAfterFix.ps1` — before/after regression wrapper
  around the same comparison; see [Regression testing](#regression-testing).
- `tools\Run-MpiTaskGridChecks.ps1` — compares single-process output against
  MPI output for the same case; see [Regression Testing
  MPI](#regression-testing-mpi).

There is currently no bundled plotting or format-conversion tool; the
tab-separated `.dat` files are read directly with any tool that accepts a
whitespace-delimited table with a header row (grid files) or without one
(particle files) — e.g. `numpy.loadtxt`, `pandas.read_csv(sep="\s+")`, or a
spreadsheet's "delimited text" import.

## 12. Recommended Simulation Cases

For the electron-thermal-transport study, use a small controlled matrix:

```text
Case A: cold electrons, stationary background
Case B: finite electron temperature, chi = 0, stationary background
Case C: finite electron temperature, chi = 1e-12 h^2/tau, stationary background
Case D: finite electron temperature, chi = 1e-11 h^2/tau, stationary background
Case E: finite electron temperature, chi = 5e-11 h^2/tau, stationary background
Case F: selected finite-temperature case with flowing background
```

Case A is the numerical and physical reference. It is expected to show stronger grid imprinting in the electric field. Cases B-D test whether electron pressure and heat transport reduce the grid imprint without applying an artificial field filter.

## 13. OpenMP Performance Testing

OpenMP scaling should be measured with deterministic particle placement and minimal diagnostics. There is currently no automated scaling harness bundled with this repository; measure scaling manually by running the same task at increasing `OMP_NUM_THREADS` and comparing wall time, optionally with a build configured with `OPENPIC_PROFILE_PARTICLE_PHASES=ON` for per-phase timing (see [CMake build options](#cmake-build-options)).

On the measured 24-logical-thread workstation, the particle half-step continues to scale up to 24 threads, but whole-run wall time saturates much earlier. `OMP_NUM_THREADS=16` is the current recommended default for efficient workstation runs; `24` gives the shortest measured wall time but with poor CPU efficiency. There is no `docs/performance/` report checked into this repository yet backing that number — treat it as a rule of thumb pending a re-measurement.

## 14. MPI Distributed Mode

### Overview

OpenPIC supports distributed execution via MPI using a **replicated-grid / distributed-particles** model. Each MPI process holds a complete copy of the simulation grid. Particles are split evenly across processes. After each push-and-scatter step, density moments are summed across all processes so that every process arrives at the same grid state and can advance the fields independently.

This model is simpler than a domain-decomposed MPI: no halo exchange, no particle migration between domains. Particle storage is rank-local for tasks that create particles through the injection callbacks (all production tasks), so MPI avoids full-particle-array memory spikes. Tasks that fill `pic_particles` directly from Lua (the `04_VERIFICATION` suite) are not rank-partitioned and are restricted to single-process runs. The tradeoff is that grid memory is not reduced: each process stores the full grid.

### When to use MPI

MPI is useful when the bottleneck is particle push throughput, not memory. A run with 10 million particles and 8 MPI ranks will scatter roughly 1.25 million particles per rank per time step instead of 10 million.

OpenMP alone is typically sufficient for moderate particle counts (up to ~5 million) on a single workstation. Add MPI when the run fills several hours and the node has multiple sockets or you have access to a cluster.

### Build with MPI support

MS-MPI (Windows) or OpenMPI (Linux) must be installed before configuring CMake.

**Windows — install MS-MPI:**
```text
powershell -ExecutionPolicy Bypass -File tools\setup_msmpi.ps1
```

**Linux — install OpenMPI:**
```bash
bash tools/setup_mpi_linux.sh
```

**CMake configuration** (same command on Windows and Linux):
```text
cmake -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release
```

A build without `-DOPENPIC_ENABLE_MPI=ON` behaves exactly as before and has no MPI dependency.

### Running with MPI

Windows:
```text
mpiexec -n 4 build\bin\Release\open-pic.exe -mpi sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Linux:
```bash
mpiexec -n 4 build/bin/open-pic -mpi sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

Without `-mpi` the single-process path is used regardless of how the binary was built.

Windows:
```text
build\bin\Release\open-pic.exe sim\tasks\04_VERIFICATION\UNIFORM_EQUILIBRIUM\main.lua
```

Linux:
```bash
build/bin/open-pic sim/tasks/04_VERIFICATION/UNIFORM_EQUILIBRIUM/main.lua
```

The number of MPI ranks (`-n`) should not exceed the number of physical CPU cores. Each rank also spawns OpenMP threads, so the total thread count is `ranks × OMP_NUM_THREADS`. Set `OMP_NUM_THREADS` explicitly to avoid oversubscription.

Windows:
```text
set OMP_NUM_THREADS=4
mpiexec -n 4 open-pic.exe -mpi main.lua   & rem 16 threads total on a 16-core machine
```

Linux:
```bash
OMP_NUM_THREADS=4 mpiexec -n 4 build/bin/open-pic -mpi main.lua   # 16 threads total on a 16-core machine
```

### What changes in MPI mode

| Feature | Single-process | MPI |
|---------|---------------|-----|
| Grid | one full copy | full copy on every rank |
| Particles | all on one process | split evenly by index |
| Field computation | all steps | identical on every rank (after allreduce) |
| `save_grid` | every save step | rank 0 only, after density reduction |
| `save_particles` | every save step | **skipped** (not yet implemented) |
| Particle diagnostics | included | **skipped** |
| Background injection | all steps | rank 0 only |
| Lua `proc_idx` | always 0 | MPI rank index |
| Lua `process_num` | always 1 | total rank count |

### How the density reduction works

After each particle push, every rank scatters its particle subset into a local density buffer. Then `MPI_Allreduce` sums the buffers across all ranks — each rank receives the global sum. The fields (B, E, UE, Te) are then computed independently on each rank from this identical density. Because the input is the same, the output is the same, so no field broadcast is needed.

The reduction covers three quantities per grid cell: `NP`, `UP` (ion bulk velocity × density), and `UP_NP` (face-density accumulator used for velocity normalization). All other cell fields are not communicated.

### Regression Testing MPI

MPI regression checks compare saved `*grd*.dat` files from a single-process run and an MPI run. This tool is PowerShell-only (Windows):

```text
powershell -ExecutionPolicy Bypass -File tools\Run-MpiTaskGridChecks.ps1
```

The test harness creates temporary short-run copies of task files and does not modify the original tasks. Passed test directories are deleted automatically; failed test directories are kept for inspection. Use `-KeepPassedArtifacts` when you need to inspect successful outputs. The harness intentionally disables `cloud_jitter_enabled` and `backgr_jitter_enabled` in those copies. Random jitter changes the particle initialization stream when particles are distributed across ranks, so it is not a stable oracle for strict serial-vs-MPI grid-file equality. Jitter physics can still be tested separately, but MPI grid equivalence tests should use deterministic particle placement.

For tests that check physical correctness rather than implementation
equivalence — analytic solutions, conservation checks, and benchmark
models — see the [Verification tasks](#verification-tasks) suite above;
there is no separate verification/validation plan document in this
repository.

### Lua scripts in MPI mode

Existing Lua scripts work without modification. The `print_mpi.print_root` calls already suppress output on non-root ranks. The `pic_parameters.process_idx` field available in Lua is set to the MPI rank.

The boundary-condition callbacks (`on_set_boundary_*`) are called on every rank because the grid is replicated. Since all ranks share the same density after reduction and compute fields identically, the boundary values are the same on every rank.

### Limitations of the current MPI mode

- **Particle diagnostics are not saved** — `save_particles` is disabled. Grid diagnostics (`save_grid`) are complete.
- **PSTD solver** — works in MPI mode (each rank runs the full FFT), but doubles compute work per rank relative to FDTD. Use FDTD for MPI production runs.
- **Grid memory is not reduced** — each rank holds the full grid. Particle memory is distributed by rank.
- **Group particle diagnostics** — cloud and background group grids are correctly reduced and saved by rank 0.

## 15. Troubleshooting

### `BoundaryKind` is nil in Lua

The executable is older than the Lua script. Rebuild the Release executable and make sure the run script uses the freshly built binary:

```text
build\bin\Release\open-pic.exe   (Windows)
build/bin/open-pic               (Linux)
```

### `opic_trace.log` grows rapidly

This usually means a diagnostic trace is being written for every low-density cell. Check whether `backgr_fracture` or another threshold function is logging per-cell messages. Production runs should not write per-cell trace messages.

### Edge density or velocity is wrong

Check the boundary-condition order and the quantity layout. `NP` is cell-centered, while `UPx`, `UPy`, and `UPz` are face-centered. The denominator used for velocity normalization must match the face location.

### Large `div(E)`

The hybrid model computes `E` from Ohm's law rather than by solving a Poisson equation. Therefore the diagnostic `div(E)` is not constrained to vanish. Large or square-shaped structures in `E` or `div(E)` should be compared across cold-electron and finite-temperature cases.

### PSTD time step is very small

PSTD has a different stability condition from FDTD. Use PSTD as a controlled numerical experiment, and keep FDTD as the reference unless the simulation series is specifically designed to compare field solvers.
