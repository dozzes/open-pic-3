# sim/tasks — Simulation Cases

Cases are grouped into subdirectories by physical scenario.  Each case
contains a single `main.lua` that sets physical and numerical parameters
and calls `require("sim_core").run()`.  Shared modules live in `../lib/`.

```
sim/tasks/
├── 01_BASELINE/    solver and lattice comparison (Ma=2, uniform B)
├── 02_BG_FLOW/     background flow variation (uniform B, high Ma)
├── 03_DIPOLE/      dipole field geometry
└── SOLAR_WIND/     background-only calibration (no cloud)
```

Run any case from its own directory:

```
cd sim/tasks/01_BASELINE/MA2_D44_UNIF_COLD
opic main.lua
```

---

## Directory Naming Convention

```
MA{N}_D{delta}_{FIELD}_{COLD|THERM}[_BG{M}][_modifiers]
```

All tags uppercase, separated by `_`.  Decimal points written as `P`.

| Tag | Meaning | Examples |
|-----|---------|---------|
| `MA{N}` | Alfven Mach number of cloud expansion | `MA2`, `MA15` |
| `D{delta}` | Magnetolaminar parameter Delta = (Rg/RL)^2 | `D44`, `D0P8`, `D7` |
| `{FIELD}` | Magnetic field geometry | `UNIF`, `DIP_POLE`, `DIP_EQUAT`, `PLANET_DIP` |
| `COLD` / `THERM` | Electron model | `COLD` = ideal (no resistivity), `THERM` = Spitzer |
| `BG{M}` | Background Alfven Mach (omit if stationary) | `BG0P5`, `BG7P3`, `BG14P5` |
| `N{exp}` | Non-default background density | `N13` = 1e13, `N14` = 1e14 cm^-3 |
| `PSTD` | Spectral field solver (default is FDTD) | |
| `JITTER` | Random offset on particle lattice | |
| `LOWRES` | Reduced resolution, quick test | |
| `X2T` | Extended duration (2x break_times) | |

Delta regimes: Delta >> 1 — electromagnetic braking; Delta << 1 — gas-dynamic braking.

---

## Groups

### 01_BASELINE — Solver and Lattice Comparison (Ma=2, uniform B)

Reference regime: low Mach number (Ma=2), high Delta (44) — deep
electromagnetic braking.  The cloud barely penetrates the background field.
These cases fix the physics and vary only the numerical algorithm, making
them the primary tool for verifying numerical convergence.

| Directory | What differs |
|-----------|-------------|
| `MA2_D44_UNIF_COLD` | Reference: FDTD solver, regular particle lattice |
| `MA2_D44_UNIF_COLD_JITTER` | Jitter ON: random offset breaks lattice periodicity |
| `MA2_D44_UNIF_COLD_JITTER_2` | Repeat of jitter run (different random seed) |
| `MA2_D44_UNIF_COLD_PSTD` | PSTD spectral solver instead of FDTD |
| `MA2_D44_UNIF_COLD_LOWRES` | Halved resolution: quick sanity check |

---

### 02_BG_FLOW — Background Flow Variation (uniform B)

Gas-dynamic regime: high Mach number (Ma=15), low Delta (~0.8) — cloud
expands faster than the field can brake it.  The background plasma flows
along X.  These cases trace how the braking and shock structure change as
the solar-wind speed increases from near-stationary to super-Alfvenic.

| Directory | Background Mach | Notes |
|-----------|----------------|-------|
| `MA8_D7_UNIF_COLD_N13_BG0P5` | 0.5 V_A | Lower cloud Mach (Ma=8), lower density (1e13) |
| `MA15_D0P8_UNIF_COLD_BG0P5` | 0.5 V_A | Slow background, near-stationary |
| `MA15_D0P8_UNIF_COLD_BG7P3` | 7.3 V_A | Intermediate background flow |
| `MA15_D0P8_UNIF_COLD_BG14P5_X2T` | 14.5 V_A | Super-Alfvenic background; 2x run duration |

---

### 03_DIPOLE — Dipole Field Geometry

Same gas-dynamic regime (Ma=15), but the background field is not uniform.
A magnetic dipole is placed outside the simulation box; the cloud expands
into the dipole fringe field.  These cases study how field-line geometry
affects cloud morphology and electromagnetic energy transfer.

| Directory | Field configuration |
|-----------|-------------------|
| `MA15_D0P8_DIP_POLE_COLD` | Cloud above magnetic pole; n=1e15 |
| `MA15_D0P7_DIP_POLE_COLD_N14` | Same geometry; lower density (1e14), Delta=0.7 |
| `MA15_D0P8_DIP_EQUAT_COLD_BG7P3` | Cloud in equatorial plane; background at 7.3 V_A |
| `MA15_D0P7_PLANET_DIP_COLD_N14` | Full planetary dipole with absorptive body; n=1e14 |

---

### SOLAR_WIND — Background-Only Calibration

**No cloud.**  A uniform plasma slab fills the box at t=0 and streams
continuously through a uniform Bz field at V_bg = 2e7 cm/s (~ 9.2 V_A).
Background density is 1e14 cm^-3; V_A0 ~ 2.18e6 cm/s.

Purpose: calibrate boundary conditions and the background inflow mechanism
in isolation, before combining with a cloud.  This is the solar-wind
streaming problem without any obstacle.

> **Status: requires extension.**  `sim_core` must support
> `background_only = true` to skip cloud initialization and fill the full
> box with background particles.  The fixed 81-node grid of the original
> run also needs an `x_nodes_override` hook.  Do not run as-is.
