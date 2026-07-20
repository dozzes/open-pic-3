# OpenPIC User Guide

## 1. Overview

OpenPIC is a hybrid particle-in-cell code for modeling ion kinetic dynamics in a magnetized plasma. Ions are represented by macroparticles, while electrons are treated as a massless fluid. The electromagnetic fields are advanced on a staggered C-grid, and the electric field is computed from a generalized Ohm's law.

The code is driven by Lua input files. A typical simulation directory contains the main physical setup, boundary-condition scripts, diagnostic settings, and a run script. The executable reads the Lua configuration, initializes particle groups and fields, advances the simulation, and writes grid and particle diagnostics.

## 2. Build

The current development setup uses CMake and Microsoft Visual Studio on Windows.

Recommended build:

```powershell
cmake --build C:\Work\open-pic-2\build --config Release
```

The Release executable is produced at:

```text
C:\Work\open-pic-2\build\bin\Release\open-pic.exe
```

For production runs, keep the compiler and build configuration fixed across a simulation series. Changing the compiler, optimization level, or OpenMP runtime can slightly change performance and numerical noise.

### CMake build options

All optional OpenPIC CMake switches are disabled by default. A plain Release
build is the production build:

```powershell
cmake -S C:\Work\open-pic-2 -B C:\Work\open-pic-2\build
cmake --build C:\Work\open-pic-2\build --config Release
```

| Option | Default | Use for | Production impact |
|---|---:|---|---|
| `OPENPIC_BUILD_TESTS` | `OFF` | Builds `openpic_tests` and enables CTest verification tests. | No effect on `open-pic.exe`; adds a separate test executable. |
| `OPENPIC_ENABLE_MPI` | `OFF` | Builds MPI-capable OpenPIC. Required for `open-pic.exe -mpi ...`. | Adds MPI dependency and MPI code paths. Without `-mpi`, the executable still uses the single-process path. |
| `OPENPIC_PROFILE_PARTICLE_PHASES` | `OFF` | Prints detailed timing blocks for particle phases and high-level simulation phases. | Profiling only. Do not use for production timing except controlled performance tests. |
| `OPENPIC_PARTICLE_STATIC_SCHEDULE` | `OFF` | Replaces the particle half-step OpenMP loop schedule with `schedule(static)` for scaling experiments. | Can change performance characteristics. Treat as experimental until validated for production workloads. |

Recommended production build:

```powershell
cmake -S C:\Work\open-pic-2 -B C:\Work\open-pic-2\build
cmake --build C:\Work\open-pic-2\build --config Release
```

MPI build:

```powershell
cmake -S C:\Work\open-pic-2 -B C:\Work\open-pic-2\build-mpi -DOPENPIC_ENABLE_MPI=ON
cmake --build C:\Work\open-pic-2\build-mpi --config Release
```

OpenMP profiling build:

```powershell
cmake -S C:\Work\open-pic-2 -B C:\Work\open-pic-2\build-omp-profile `
  -DOPENPIC_PROFILE_PARTICLE_PHASES=ON `
  -DOPENPIC_PARTICLE_STATIC_SCHEDULE=ON `
  -DOPENPIC_ENABLE_MPI=OFF
cmake --build C:\Work\open-pic-2\build-omp-profile --config Release
```

Do not mix profiling builds with production baselines. The profiling code is
compiled out unless `OPENPIC_PROFILE_PARTICLE_PHASES=ON`, but when enabled it
adds timers and log output. `OPENPIC_PARTICLE_STATIC_SCHEDULE=ON` is also a
deliberate performance experiment, not the default production path.

### Runtime OpenMP setting

`OMP_NUM_THREADS` is not a CMake option. It is a runtime environment variable
that controls how many OpenMP worker threads the executable uses:

```powershell
$env:OMP_NUM_THREADS = 16
C:\Work\open-pic-2\build\bin\Release\open-pic.exe main.lua
```

On the measured 24-logical-thread workstation, `16` is the recommended practical
setting for single-process OpenMP runs. `24` can be slightly faster in wall time
but uses CPU less efficiently.

### Build and run C++ verification tests

OpenPIC has a lightweight C++ verification test executable that does not depend
on GoogleTest or other external test frameworks. It is disabled by default and
can be enabled with:

```powershell
cmake -S C:\Work\open-pic-2 -B C:\Work\open-pic-2\build-tests -DOPENPIC_BUILD_TESTS=ON
cmake --build C:\Work\open-pic-2\build-tests --config Release --target openpic_tests
ctest --test-dir C:\Work\open-pic-2\build-tests -C Release --output-on-failure
```

The first test target is `openpic_tests`. It currently checks:

- CIC scatter weight conservation for `NP`, `UP`, and `UP_NP`.
- NGP and quadratic TSC gather/deposit consistency on the staggered grid.
- Boris pusher in uniform electric field against the analytic velocity update.
- Boris pusher in uniform magnetic field for speed conservation.

### Run verification benchmark tasks

The Lua-level verification and benchmark smoke suite is run with:

```powershell
powershell -ExecutionPolicy Bypass -File tools\Run-VerificationBenchmarks.ps1
```

Run a subset with:

```powershell
powershell -ExecutionPolicy Bypass -File tools\Run-VerificationBenchmarks.ps1 `
    -Case uniform-equilibrium,linear-alfven,field-loop-divb
```

For the list of implemented cases, what each one checks, and the operational
details (`smoke_ok` vs `pass`, `-KeepArtifacts`, run-directory cleanup), see
`docs/Verification_Validation_Plan.md`.

## 3. Run Directory Layout

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

## 4. Lua Configuration

The main simulation parameters are defined in `main.lua`.

Common physical parameters:

```lua
Bz0          -- initial magnetic field
backgr_dens  -- background plasma density
cloud_ions_num
cloud_R
cloud_v
```

Common numerical parameters:

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
openpic_tests: test_pstd_plane_wave_faraday_derivative
tools/Run-VerificationBenchmarks.ps1 -Case pstd-baseline-smoke
```

This unit test initializes a sinusoidal `Ez(y + h/2)` mode and checks that the
PSTD Faraday half-step updates `Bx` with the analytic spectral derivative.
The Lua smoke case loads the real `MA2_D44_UNIF_COLD_PSTD` task but caps it at
one step. A two-step run of that baseline currently fails with a particle CFL
violation after spectral-field runaway, so do not treat the PSTD baseline as a
validated production case yet.

## 5. Electron Model

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

## 6. Particle Groups

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

## 7. Boundary Conditions

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

## 8. Stationary and Flowing Backgrounds

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

## 9. Diagnostics

Grid diagnostics are tab-separated files with columns:

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

## 10. Recommended Simulation Cases

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

## 11. OpenMP Performance Testing

OpenMP scaling should be measured with deterministic particle placement and minimal diagnostics. The current profiling harness is:

```powershell
powershell -ExecutionPolicy Bypass -File tools\Run-OpenMpScaling.ps1 `
  -TaskMain sim\tasks\01_BASELINE\MA2_D44_UNIF_COLD\main.lua `
  -Threads 1,2,4,8,12,16,24 `
  -TimeSteps 10 `
  -CloudPartsOnStep 25 `
  -BackgrPartsOnStep 4 `
  -KeepArtifacts
```

The harness writes `openmp_scaling.csv` and parses both high-level simulation timings and particle-loop phase timings when the binary is built with `OPENPIC_PROFILE_PARTICLE_PHASES=ON`.

The current saturation report is in `docs/performance/OpenMP_Performance_Report.md`. On the measured 24-logical-thread workstation, the particle half-step continues to scale up to 24 threads, but whole-run wall time saturates much earlier. `OMP_NUM_THREADS=16` is the current recommended default for efficient workstation runs; `24` gives the shortest measured wall time but with poor CPU efficiency.

## 12. MPI Distributed Mode

### Overview

OpenPIC supports distributed execution via MPI using a **replicated-grid / distributed-particles** model. Each MPI process holds a complete copy of the simulation grid. Particles are split evenly across processes. After each push-and-scatter step, density moments are summed across all processes so that every process arrives at the same grid state and can advance the fields independently.

This model is simpler than a domain-decomposed MPI: no halo exchange, no particle migration between domains. Particle storage is rank-local for tasks that create particles through the injection callbacks (all production tasks), so MPI avoids full-particle-array memory spikes. Tasks that fill `pic_particles` directly from Lua (the `04_VERIFICATION` suite) are not rank-partitioned and are restricted to single-process runs. The tradeoff is that grid memory is not reduced: each process stores the full grid.

### When to use MPI

MPI is useful when the bottleneck is particle push throughput, not memory. A run with 10 million particles and 8 MPI ranks will scatter roughly 1.25 million particles per rank per time step instead of 10 million.

OpenMP alone is typically sufficient for moderate particle counts (up to ~5 million) on a single workstation. Add MPI when the run fills several hours and the node has multiple sockets or you have access to a cluster.

### Build with MPI support

MS-MPI (Windows) or OpenMPI (Linux) must be installed before configuring CMake.

**Windows — install MS-MPI:**
```powershell
powershell -ExecutionPolicy Bypass -File tools\setup_msmpi.ps1
```

**Linux — install OpenMPI:**
```bash
bash tools/setup_mpi_linux.sh
```

**CMake configuration:**
```powershell
cmake -B build -DOPENPIC_ENABLE_MPI=ON
cmake --build build --config Release
```

A build without `-DOPENPIC_ENABLE_MPI=ON` behaves exactly as before and has no MPI dependency.

### Running with MPI

```powershell
mpiexec -n 4 build\bin\Release\open-pic.exe -mpi sim\tasks\01_BASELINE\MA2_D44_UNIF_COLD\main.lua
```

Without `-mpi` the single-process path is used regardless of how the binary was built:
```powershell
build\bin\Release\open-pic.exe sim\tasks\01_BASELINE\MA2_D44_UNIF_COLD\main.lua
```

The number of MPI ranks (`-n`) should not exceed the number of physical CPU cores. Each rank also spawns OpenMP threads, so the total thread count is `ranks × OMP_NUM_THREADS`. Set `OMP_NUM_THREADS` explicitly to avoid oversubscription:
```powershell
$env:OMP_NUM_THREADS = 4
mpiexec -n 4 open-pic.exe -mpi main.lua   # 16 threads total on a 16-core machine
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

MPI regression checks compare saved `*grd*.dat` files from a single-process run and an MPI run. Use:

```powershell
powershell -ExecutionPolicy Bypass -File tools\Run-MpiTaskGridChecks.ps1
```

The test harness creates temporary short-run copies of task files and does not modify the original tasks. Passed test directories are deleted automatically; failed test directories are kept for inspection. Use `-KeepPassedArtifacts` when you need to inspect successful outputs. The harness intentionally disables `cloud_jitter_enabled` and `backgr_jitter_enabled` in those copies. Random jitter changes the particle initialization stream when particles are distributed across ranks, so it is not a stable oracle for strict serial-vs-MPI grid-file equality. Jitter physics can still be tested separately, but MPI grid equivalence tests should use deterministic particle placement.

For tests that check physical correctness rather than implementation
equivalence, see `docs/Verification_Validation_Plan.md`. Those tests should
cover analytic solutions, conservation checks, and benchmark models.

### Lua scripts in MPI mode

Existing Lua scripts work without modification. The `print_mpi.print_root` calls already suppress output on non-root ranks. The `pic_parameters.process_idx` field available in Lua is set to the MPI rank.

The boundary-condition callbacks (`on_set_boundary_*`) are called on every rank because the grid is replicated. Since all ranks share the same density after reduction and compute fields identically, the boundary values are the same on every rank.

### Limitations of the current MPI mode

- **Particle diagnostics are not saved** — `save_particles` is disabled. Grid diagnostics (`save_grid`) are complete.
- **PSTD solver** — works in MPI mode (each rank runs the full FFT), but doubles compute work per rank relative to FDTD. Use FDTD for MPI production runs.
- **Grid memory is not reduced** — each rank holds the full grid. Particle memory is distributed by rank.
- **Group particle diagnostics** — cloud and background group grids are correctly reduced and saved by rank 0.

## 13. Troubleshooting

### `BoundaryKind` is nil in Lua

The executable is older than the Lua script. Rebuild the Release executable and make sure the run script uses:

```text
C:\Work\open-pic-2\build\bin\Release\open-pic.exe
```

### `opic_trace.log` grows rapidly

This usually means a diagnostic trace is being written for every low-density cell. Check whether `backgr_fracture` or another threshold function is logging per-cell messages. Production runs should not write per-cell trace messages.

### Edge density or velocity is wrong

Check the boundary-condition order and the quantity layout. `NP` is cell-centered, while `UPx`, `UPy`, and `UPz` are face-centered. The denominator used for velocity normalization must match the face location.

### Large `div(E)`

The hybrid model computes `E` from Ohm's law rather than by solving a Poisson equation. Therefore the diagnostic `div(E)` is not constrained to vanish. Large or square-shaped structures in `E` or `div(E)` should be compared across cold-electron and finite-temperature cases.

### PSTD time step is very small

PSTD has a different stability condition from FDTD. Use PSTD as a controlled numerical experiment, and keep FDTD as the reference unless the simulation series is specifically designed to compare field solvers.
