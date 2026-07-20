# sim/lib — OpenPIC 3D Simulation Bundle

## Overview

This directory is a **reusable simulation bundle** for the OpenPIC 3D hybrid
particle-in-cell (PIC) code.  It implements a plasma cloud expanding into a
magnetized background plasma — a canonical problem in space and laboratory
plasma physics (comets, laser-ablated pellets, artificial plasma releases).

The bundle separates concerns strictly:

- **One file per concern** — physics, particles, fields, diagnostics, logging
  are each in their own module.
- **One file per case** — to define a new simulation, write only a `main.lua`
  that sets scenario parameters and calls `require("sim_core").run()`.
- **All shared modules live here** — cases in sibling directories point their
  `package.path` at this directory and inherit everything automatically.

---

## Directory Layout

```
sim/
|-- lib/                    This directory (shared modules).
|   |-- README.md           This file.
|   |-- cgs.lua             Physical constants in CGS-Gaussian units.
|   |-- sim_core.lua        Initialization orchestrator.  Entry point called by main.lua.
|   |-- physics.lua         Plasma physics derivation (frequencies, scales, CFL, grid).
|   |-- init_particles.lua  Particle placement: cloud sphere + background lattice.
|   |-- init_field.lua      EM field and electron temperature/resistivity initialization.
|   |-- callbacks.lua       Engine callback functions (boundary conditions, diagnostics).
|   |-- diag.lua            Per-step energy logging and .dat file archiving.
|   |-- params_log.lua      Full parameter dump to lua_params.txt.
|   |-- bound_cond.lua      Boundary condition implementations (staggered C-grid aware).
|   |-- em_field.lua        EM field initialization utilities (uniform, dipole, ...).
|   |-- parts_count.lua     Particle count estimation (cloud sphere, background box).
|   |-- print_mpi.lua       MPI-aware console output (prints only from rank 0).
|   |-- utils.lua           Portable OS utilities (mkdir, copy, remove, clear screen).
|
|-- tasks/                  One subdirectory per simulation case.
    |-- MA2_D44_UNIF_COLD/
    |-- MA2_D44_UNIF_COLD_JITTER/
    |-- MA15_D0P8_DIP_POLE_COLD/
    |-- ...
```

---

## Case Naming Convention

Cases live in `sim/tasks/`.  Each subdirectory name encodes the key physical
parameters so that the run is identifiable without opening `main.lua`:

```
MA{N}_D{delta}_{field}_{cold|therm}[_BG{M}][_modifiers]
```

All tags are **UPPERCASE**, separated by `_`.  Decimal points are written as
`P` (e.g. `D0P8` = Delta 0.8, `BG7P3` = background Mach 7.3).

| Tag | Meaning | Examples |
|-----|---------|---------|
| `MA{N}` | Alfven Mach number of cloud | `MA2`, `MA15` |
| `D{delta}` | Magnetolaminar parameter Delta = (Rg/RL)^2 | `D44`, `D0P8`, `D7` |
| `{field}` | Magnetic field geometry | `UNIF`, `DIP_POLE`, `DIP_EQUAT`, `PLANET_DIP` |
| `COLD` / `THERM` | Electron model | `COLD` = ideal (no resistivity), `THERM` = Spitzer |
| `BG{M}` | Background Alfven Mach (omit if stationary) | `BG0P5`, `BG7P3`, `BG14P5` |
| `N{exp}` | Non-default background density | `N13` = 1e13 cm^-3, `N14` = 1e14 cm^-3 |
| `PSTD` | Spectral (PSTD) field solver (default is FDTD) | `PSTD` |
| `JITTER` | Random offset on particle lattice | `JITTER` |
| `LOWRES` | Reduced resolution test run | `LOWRES` |
| `X2T` | Extended duration (2x `break_times`) | `X2T` |

**Delta regimes:**
- Delta >> 1 — electromagnetic (magnetolaminar) braking
- Delta << 1 — gas-dynamic braking
- Delta = (Rg/RL)^2 where Rg is the magnetic braking radius, RL is the Larmor radius

---

## Quick Start: Creating a New Case

1. Create a directory under `sim/tasks/`, named by the convention above.
2. Write `sim/tasks/MY_CASE/main.lua` with the structure below.
3. Run: `opic main.lua` from inside `sim/tasks/MY_CASE/`.

```lua
-- sim/tasks/MY_CASE/main.lua
package.path = package.path .. ";../../lib/?.lua"
local cgs = require("cgs")

case_name = "My case: description"

-- Plasma parameters
backgr_mi       = 1 * cgs.mp
backgr_dens     = 1.0e+15
backgr_Z        = 1
Bz0             = 100.0
cloud_mi        = 1 * cgs.mp
cloud_ions_num  = 6.0e+18
cloud_Z         = 1
Ma              = 2.0
backgr_flow_Ma_x = 0.0
backgr_flow_Ma_y = 0.0
backgr_flow_Ma_z = 0.0

-- Algorithm
magnetic_field_alg     = MagneticFieldAlg.FDTD
cold_electrons_enabled = true
cloud_jitter_enabled   = false
backgr_jitter_enabled  = false
pstd_quick_test_steps  = nil

-- Resolution
points_on_RL         = 10
h_R0                 = 4.0
break_times          = 4.0
cloud_parts_on_step  = 25
backgr_parts_on_step = 4

require("sim_core").run()
```

---

## Module Reference

### cgs.lua

Defines physical constants in CGS-Gaussian units as a plain table.

**Why CGS-Gaussian?**
Plasma physics historically uses CGS-Gaussian because Maxwell's equations take
a symmetric form and the speed of light `c` appears explicitly in the field
coupling, making electromagnetic normalization transparent.  The relation
`E [statV/cm]`, `B [Gauss]`, `n [cm^-3]` is self-consistent without
conversion factors.

| Symbol | Value | Units | Meaning |
|--------|-------|-------|---------|
| `c`    | 2.9979e+10 | cm/s | Speed of light |
| `e`    | 4.8032e-10 | esu  | Elementary charge |
| `me`   | 9.1e-28    | g    | Electron mass |
| `mp`   | 1.6726e-24 | g    | Proton mass |

---

### sim_core.lua

**Role:** Orchestrator.  Called once by `main.lua` at startup.

Executes initialization in the following fixed order:

| Step | Call | Purpose |
|------|------|---------|
| 1 | `require("callbacks")` | Register engine callbacks before any computation |
| 2 | `require("physics").compute()` | Derive all plasma scales, CFL, grid dimensions |
| 3 | `require("params_log").write(f)` | Write lua_params.txt for reproducibility |
| 4 | Memory check | Warn user if estimated RAM exceeds 2 GB |
| 5 | `require("init_particles").init()` | Place cloud and background macro-particles |
| 6 | `require("init_field").init()` | Set B, E, Te, eta on the grid |
| 7 | `pic_parameters.*` assignment | Configure solver algorithms and output scales |

The ordering is not arbitrary.  In particular, callbacks must be registered
before step 2 because the engine may trigger boundary conditions as soon as the
grid is resized; and the grid must be sized (step 2) before any particle
placement (step 5).

---

### physics.lua

**Role:** Derives every quantity that depends on the scenario parameters.
Takes the raw user inputs from `main.lua` and computes the full set of plasma
physics variables, numerical parameters, and grid geometry.

**Key quantities derived:**

| Variable | Formula | Physical meaning |
|----------|---------|-----------------|
| `V_A0` | `B / sqrt(4*pi*n*m)` | Alfven speed — the speed of low-frequency EM waves |
| `backgr_Wpi` | `Z*e * 2*sqrt(pi*n/m)` | Ion plasma frequency |
| `backgr_c_Wpi` | `c / Wpi` | Ion inertial length (skin depth) |
| `backgr_Wci` | `Z*e*B / (m*c)` | Ion cyclotron frequency |
| `RL` | `m*v*c / (Z*e*B)` | Ion Larmor radius |
| `h` | `max(2*c/Wpi, RL/N)` | Grid step — resolves both inertial length and Larmor orbit |
| `tau` | `0.2 * min(CFL conditions)` | Time step with 5x safety margin |
| `Rb` | `(6*W0 / B^2)^(1/3)` | Magnetic braking radius |
| `Tt` | `min(Tb, Tg)` | Actual braking time (sets simulation duration) |

**CFL conditions** (all must be satisfied simultaneously):

| Condition | Formula | Constraint |
|-----------|---------|-----------|
| `CFL_tau` | `0.5 * h / V_max` | Wave must not cross more than one cell per step |
| `PSTD_tau` | `0.367 * h / V_A0` | Alfven CFL for the spectral Faraday solver |
| `Wci_tau` | `0.1 / Wci_max` | At least 10 time steps per cyclotron period |
| `CFL_tau_loc` | `(h/c_Wpi)^2 / (pi*sqrt(3)*Wci)` | Suppresses dispersive whistler instability |
| `CFL_tau_Pe` | `0.5 * h / V_Pe` | Electron thermal speed (hot-electron runs only) |

All derived variables are written to the global Lua environment (`_G`) so that
`bound_cond.lua`, `diag.lua`, and `params_log.lua` can access them without
explicit argument passing.

---

### init_particles.lua

**Role:** Places all macro-particles at their initial positions and velocities.

**What is a macro-particle?**
We cannot track every real ion (a cloud of 10^18 particles).  Instead, each
macro-particle represents `ni = N_real / N_macro` real ions and carries the
same charge-to-mass ratio, so it moves identically under the Lorentz force.
The statistical noise in density and momentum scales as `1/sqrt(N_macro)` per
cell — more macro-particles per cell gives a smoother result.

**Cloud initialization:**
Macro-particles sit on a regular 3D lattice with spacing `h / cloud_parts_on_step`.
The code exploits **8-octant symmetry**: only the `x >= 0, y >= 0, z >= 0`
octant is iterated explicitly; the remaining 7 octants are filled by mirroring
positions and flipping the corresponding velocity components.  This gives exact
8-fold symmetry at `t = 0` and reduces loop count by 8×.

Each cloud particle's velocity is radially outward and proportional to
displacement from the center (Hubble-law expansion):
```
v = V_max * (x, y, z) / R0
```

**Background initialization:**
Fills the domain uniformly on the same lattice, excluding:
- Positions inside the cloud sphere (already occupied by cloud particles)
- A thin boundary strip (`backgr_particle_boundary_width_*`) where particles
  would be immediately absorbed by the wall

**Particle weighting (NGP, CIC, or TSC):**

`sim_core.run({ scatter_alg = ScatterAlg.NGP })` selects nearest-grid-point
weighting consistently for both particle-to-grid deposition and
grid-to-particle field gathering. The default is `ScatterAlg.Standard`, i.e.
CIC, for backward compatibility.

`sim_core.run({ scatter_alg = ScatterAlg.TSC })` selects the quadratic
triangular-shaped-cloud B-spline. TSC uses three points per dimension on each
quantity's own staggered Yee lattice. Deposit and gather use the same shape,
so their adjoint relationship is retained away from physical boundaries.

**Controlled sub-cell translation:**

`sim_core.run({ cloud_center_shift_h = {x=0.25, y=0.25, z=0.0} })` translates
the complete initial particle configuration by the given fractions of `h`.
The cloud radius, density, velocity profile, particle weights, and counts are
unchanged. This option is intended for grid-locking diagnostics.

**CIC weighting (Cloud-in-Cell):**
During the simulation, each macro-particle scatters its charge and momentum to
the 8 surrounding grid cells using trilinear (CIC) weighting.  The weight for
cell `(i,j,k)` is proportional to the volume of the cell diagonally opposite
to the particle position.  This is the inverse of the gather operation used to
interpolate fields back to particle positions.

---

### init_field.lua

**Role:** Sets the initial electromagnetic field and electron state on the grid.

**Fields initialized:**

| Field | Value | Location on Yee grid |
|-------|-------|----------------------|
| `B` | `(0, 0, Bz0)` uniform | Cell edges (staggered) |
| `E` | `(Ex0, Ey0, 0)` motional | Cell faces (staggered) |
| `Te` | `cloud_Te` inside R0, `backgr_Te` outside | Cell center |
| `eta` | `resistivity` (0 for cold electrons) | Cell center |

**Staggered Yee C-grid:**
Fields are not co-located — each component lives at a different point in the
cell.  This staggering is the key to the FDTD method's exact discrete curl
operators and absence of spurious monopoles.  `em_field.uniform()` handles the
staggering internally.

**Motional electric field:**
When the background plasma drifts with velocity `v`, the frozen-in field
condition `E = -v x B / c` requires a non-zero electric field even in the
initial uniform state.  For flow along X in field along Z:
`Ex = -vy*Bz/c`, `Ey = vx*Bz/c`.

---

### callbacks.lua

**Role:** Defines all global functions that the C++ engine calls during the
simulation loop.

**Callback execution order per time step:**

```
on_iteration_begin()           -- inject inflow particles (if background flows)
  |
  +-- [engine] update B (half step, Faraday)
  +-- [engine] compute UE from curl(B)
  +-- [engine] compute E from Ohm's law
  +-- [engine] push particles (Boris)
  |
on_particles_moved_half_time()
  |
  +-- [engine] scatter density and momentum to grid
  |
on_set_boundary_NP()           -- fix ion density at domain walls
on_set_boundary_group_NP()     -- per-group variant
on_set_boundary_UP()           -- fix ion bulk velocity
on_set_boundary_group_UP()     -- per-group variant
on_set_boundary_UE()           -- fix electron velocity
on_set_boundary_Te()           -- fix electron temperature
on_set_boundary_eta()          -- fix resistivity
on_set_boundary_EF()           -- fix electric field
on_set_boundary_MF()           -- fix magnetic field
  |
  +-- [engine] update B (second half step, corrector)
  +-- [engine] move particles to full step
  |
on_particles_moved_full_time()
on_iteration_end()             -- energy diagnostics, file archiving
```

**Background inflow:**
For a flowing background (`backgr_flow_Ma_x != 0`), particles stream out of
the downstream boundary.  `on_iteration_begin()` injects new particle planes
at the inflow boundary to maintain a steady-state density.  The injection rate
is controlled by accumulating `vx * tau` each step and adding one plane per
`backgr_part_dist` of accumulated displacement.

---

### diag.lua

**Role:** Per-step energy diagnostics and output file management.

**Energy log (energy.txt):**

Each row corresponds to one save step and contains:

| Column | Formula | Physical meaning |
|--------|---------|-----------------|
| `t` | `step * tau / T_scale` | Normalized time (in units of 1/Wci) |
| `W_cloud` | (reserved) | Cloud ion kinetic energy / W0 |
| `W_backgr` | (reserved) | Background ion kinetic energy / W0 |
| `W_mf` | `sum(B^2) / (8*pi) / W0` | Magnetic field energy / W0 |
| `W_e` | `sum(n*me*Ue^2/2) / W0` | Electron kinetic energy / W0 |

`W0 = 0.3 * M_cloud * V_max^2` is the initial cloud kinetic energy used as
the normalization reference.

**File archiving:**
The engine writes one `.dat` file per saved field quantity per time step.  A
long run produces thousands of files.  At each save step, `diag.lua` copies
the relevant cross-section planes into `diag/`, then packs all `.dat` files
into a timestep-named `.zip` archive and removes the originals, keeping the
working directory clean.

---

### params_log.lua

**Role:** Writes a complete, annotated record of all derived simulation
parameters to `lua_params.txt` in the run directory.

Every entry includes physical units and, where applicable, a constraint that
must be satisfied for the run to be valid (e.g. `cloud_Wci*tau < 0.2`).
This file is the primary source of truth for post-processing and for
reproducing a run from scratch.

---

### bound_cond.lua

**Role:** Implements boundary condition updates for the staggered C-grid.

Provides `nonperturbed_*` functions that copy interior values into boundary
ghost cells, making the wall look like an undisturbed plasma (transparent /
open boundary).  The `nonperturbed_flow_*` variants additionally force boundary
cells to the background drift value to sustain a stationary inflow condition.

All functions are staggering-aware: face-centered components (`E`, `UP`, `UE`,
`B`) use a different extrapolation offset than cell-centered scalars (`NP`,
`Te`, `eta`).

---

### em_field.lua

**Role:** Utility library for initializing electromagnetic fields on the grid.

| Function | Description |
|----------|-------------|
| `uniform(grid, Ex,Ey,Ez,Bx,By,Bz)` | Set uniform E and B everywhere |
| `dipole_equator(grid, B0, ratio)` | Set a magnetic dipole field along the equator |

---

### parts_count.lua

**Role:** Pre-counts macro-particles before allocating the particle array.

| Function | Returns |
|----------|---------|
| `get_cloud_parts_num(R0, dist)` | Number of cloud macro-particles in the 8-octant sphere |
| `get_backgr_parts_num(...)` | Number of background macro-particles in the box minus sphere |

These functions replicate the loop structure of `init_particles.lua` without
actually placing particles, so the total count is known before `pic_particles.size`
is set.

---

### print_mpi.lua

**Role:** MPI-safe console output.

In a multi-process run each rank would otherwise print identical messages,
flooding the terminal.  `print_mpi.print_root(proc_idx, ...)` prints only when
`proc_idx == 0` (the master rank).

---

### utils.lua

**Role:** Portable OS utilities that work on both Windows and Linux/macOS.

| Symbol | Windows | Linux/macOS |
|--------|---------|-------------|
| `utils.CP` | `copy /y` | `cp` |
| `utils.RM` | `del /q` | `rm -f` |
| `utils.MOVE` | `move` | `mv` |
| `utils.clear_cmd()` | `cls` | `clear` |
| `utils.ensure_dir(d)` | `if not exist d mkdir d` | `mkdir -p d` |
| `utils.copy_pattern(p,d)` | `copy /y p d` | `cp p d` |

---

## Physics Background

### Hybrid PIC model

The code uses a **hybrid kinetic-fluid** model:

- **Ions** are kinetic: each macro-particle is pushed by the full Lorentz force
  using the Boris algorithm (time-reversible, energy-conserving).
- **Electrons** are a massless fluid: their velocity is determined
  instantaneously from `curl(B) = (4*pi/c)*J` (Ampere's law without
  displacement current, valid for `omega << Wpe`).

The field equations are:

```
Faraday:   dB/dt = -c * curl(E)
Ohm:       E = (UE x B)/c - grad(Pe)/(n*e) + eta*J
Pe = n*k_B*Te   (electron pressure, adiabatic or isothermal)
```

Ions provide `n` and `J_i` via scatter; electrons close the system.

### Time integration (leapfrog)

Fields and particles are staggered in time:

```
t:       B^(n-1/2) ... B^(n+1/2) ... B^(n+1) ...
         E^n              E^(n+1)
         r^n                 r^(n+1)
         v^(n-1/2)              v^(n+1/2)
```

The half-step stagger means the Boris push is second-order accurate in time
without requiring a corrector.

### Electron closure parameters (eta, chi, kappa)

The electron fluid is closed by three tunable coefficients. Defaults are set
inside `sim_core.run()`; case overrides go in `main.lua` AFTER the `run()`
call (see the template's "ELECTRON CLOSURE OVERRIDES" section).
Textbook-style guide with full derivations: `docs/Eta_Chi_Kappa_Guide.md`.
Detailed measurements behind every claim here: `docs/M20D14_Analysis_Notes.md`
(2026-07-08 reference entry). Example numbers below are for the
`06_ARTICLE_MA20_D14` grid (h = 1.46 cm, tau = 1.04e-8 s, n0 = 1e15 cm^-3,
B0 = 100 G, Te_ref = 1 eV).

| Parameter | Units | Acts on | Measured effect on the m=4 grid mode |
|-----------|-------|---------|--------------------------------------|
| `resistivity` (+ Spitzer knobs) | s | B (diffusion), Te (Joule) | dominant: adaptive eta(Te) -> Q4 ~ 0.018, const eta -> ~0.010 |
| `electron_thermal_conductivity` | cm^2/s | Te only | none, under any eta model |
| `kappa_friction` + `kappa_B_ref` | dimensionless / G | B, Te, ion velocities | strongly AMPLIFIES at kappa=0.1: peak Q4 2x CHI0, 5x ETACONST (measured) |

**eta — resistivity [s].**
Electron-ion collisions. Computed physically (Spitzer at `Te_ref` = 1 eV) in
`physics.lua`; do not set by hand. Enters Ohm's law as `eta*(c/4pi)*curl(B)`
and Joule heating `dTe = eta*j^2*tau/(1.5*n*k_B)`. In the induction equation
it acts as magnetic diffusion with per-step diffusion number
`D = eta*c^2*tau/(4*pi*h^2)` (~0.04 at base_eta) — the grid-scale (2h)
checkerboard in B is damped by `|1 - 12*D|` ~ 0.52 per step, making eta the
main damper of grid-scale field noise. The explicit term requires D <= 1/6;
the engine caps eta at `eta_stab` (D = 0.1) automatically. The adaptive
scaling `eta = base_eta*(Te_ref/Te)^1.5` weakens damping at the hot
compressed front (Te > Te_ref -> eta down to 0.5-0.8x base_eta) exactly
where the m=4 mode grows — the measured driver of the mode. Knobs:

- `Spitzer_eta_floor_mult = 1.0` — floor eta at base_eta; structurally binds
  only where Te > Te_ref (the front), leaves the cold cavity's adaptive
  behavior untouched. Minimal-footprint fix; suppresses the m=4 distortion
  x1.95-2.7 across all fields (NP, B, E, UP, UE).
- `Spitzer_Te_ref = 0` — constant eta everywhere (diagnostic switch).
  WARNING: the grid dump normalizes its Te column by this value, so Te
  becomes `inf` in the dump (artifact, not physics).
- `Spitzer_Te_smooth_passes` — binomial smoothing of the Te fed to eta;
  measured to have no effect (the front eta deficit is a smooth trend over
  tens of cells, not single-cell noise).
- `hyper_resistivity` (eta4 [s*cm^2]) — separate k^4 induction filter
  outside Ohm's law; selectively damps 2h modes, particles never feel it.

**chi — electron heat conduction [cm^2/s].**
`electron_thermal_conductivity`; physical Spitzer value for this regime is
`chi_e = 3.2*v_Te^2*tau_e = 1.93e5` (Te = 1 eV, n = 1e15, NRL formulary).
Solved by an implicit ADI/Thomas solver (`src/heat_solver.h`) —
unconditionally stable at any value; `<= 0` disables it entirely.
IMPLEMENTATION NOTE: the effective diffusivity is
`chi * clamp(Te/(eta + 1e-12), 1e-6, 100)`; at typical scales Te/eta ~ 1e17,
so the clamp is ALWAYS saturated — the solver effectively runs with a
uniform `100*chi`, and the advertised Spitzer-like Te^2.5 modulation never
engages. Measured repeatedly (CHI1E3, CHI_SPITZER, ETACONST_CHISPITZER):
no effect on the m=4 mode — the mode lives in B/currents, chi only diffuses
Te, and the Te->eta(Te) path is a smooth large-scale trend that per-step
diffusion (~0.09 h^2/tau even with the x100) barely changes.

**kappa — legacy-code-style interspecies friction (dimensionless).**
`kappa_friction` (analog of the legacy code's `xappa`) with `kappa_B_ref`
(analog of its field unit h0, i.e. Bz0). When both > 0 this REPLACES the
Spitzer eta model with the legacy two-fluid friction closure — three
channels as one package:

1. Ohm's law: per-cell `eta = kappa*B_ref/(c*e*ne)` — no Te feedback,
   proportional to 1/ne (weak in the dense shell, strong in rarefaction,
   capped at `eta_stab` where 1/ne diverges toward the cavity).
2. Friction heating: the same eta feeds the Joule term (equals the legacy
   `0.5*tau*xappa*|Ui-Ue|^2`).
3. Ion drag: every ion feels `du/dt = kappa*Omega_ci*(Ue - u)`, applied as
   an exact implicit sub-step after the Boris push (unconditionally stable).
   This is the ion-side half of the friction that standard hybrid schemes
   drop; it damps currents (including grid-scale current noise) directly in
   velocity space.

Conversions: `eta_equiv(n) = kappa*B_ref/(c*e*n)`;
`kappa_equiv(eta) = eta*c*e*n/B_ref`. For this regime: legacy `kappa = 0.1`
corresponds to `eta = 6.94e-16 s` at n0 (~165x weaker than Spitzer
base_eta); conversely our base_eta corresponds to `kappa ~ 16.5`
(physically `nu = kappa*Omega_ci ~ (m_e/m_i)*nu_ei`). Measured impact on
the m=4 mode (MA20_D14_UNIF_TE_KAPPA01, 2026-07-08): the damping deficit
dominates — peak Q4 = 0.061, 2x worse than CHI0 and 5x worse than
ETACONST; neither removing the Te feedback nor the ion drag at
0.1*Omega_ci compensates the ~165x weaker B-noise damping. Treat the
kappa mode as a legacy-comparison control point, not a production
setting. Keep `Spitzer_Te_ref` nonzero: the kappa branch ignores it for
eta, but the Te dump column is still normalized by it.

### Grid staggering (Yee C-grid)

```
Cell (i,j,k):

  B components (edges):
    Bx at (i+1/2, j,     k    )
    By at (i,     j+1/2, k    )
    Bz at (i,     j,     k+1/2)

  E components (faces):
    Ex at (i,     j+1/2, k+1/2)
    Ey at (i+1/2, j,     k+1/2)
    Ez at (i+1/2, j+1/2, k    )

  Scalars (cell center):
    NP, Te, eta at (i+1/2, j+1/2, k+1/2)
```

The staggering ensures that the discrete curl operators satisfy
`div(curl(E)) = 0` exactly, preventing spurious charge accumulation.
