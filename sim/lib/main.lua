--[[
  main.lua -- OpenPIC 3D simulation entry point
  =============================================
  This is the ONLY file you need to edit to set up a new simulation case.
  It defines:
    1. Physical scenario  -- what plasma we are simulating
    2. Algorithm choices  -- which numerical methods to use
    3. Resolution         -- accuracy vs. cost trade-offs

  Everything else (physics derivation, particle placement, field setup,
  diagnostics) lives in separate modules and is called at the bottom via
  require("sim_core").run().

  HOW A HYBRID PIC SIMULATION WORKS (brief overview)
  ---------------------------------------------------
  We model two populations:
    * Ions  -- treated as individual particles (kinetic)
    * Electrons -- treated as a fluid (fast, respond instantly)

  At each time step the code does:
    1. Update magnetic field B  (Faraday's law, half step)
    2. Compute electron velocity UE  from curl(B) and ion density
    3. Compute electric field E  from generalized Ohm's law
    4. Push ions with the Boris algorithm  (Lorentz force)
    5. Scatter ion density and momentum to the grid (selected CIC/NGP weighting)
    6. Update B to the full step, move ions to full step

  The scenario we simulate here: a plasma cloud (e.g. a comet or
  laser-ablated pellet) expanding into a background magnetized plasma.
  The cloud brakes against the background magnetic field and generates
  Alfven waves, bow shocks, and non-linear structures.
]]

package.path = package.path .. ";./?.lua"
local cgs = require("cgs")   -- physical constants in CGS-Gaussian units

-- ===========================================================================
-- CASE DESCRIPTION
-- ===========================================================================
-- Give every run a human-readable label so you can identify output files.
case_name = "Case A: cold electrons, stationary background, FDTD solver"

-- ===========================================================================
-- BACKGROUND PLASMA
-- ===========================================================================
-- The background is a uniform, magnetized, stationary (or slowly drifting)
-- hydrogen plasma that the cloud expands into.

backgr_mi   = 1 * cgs.mp    -- background ion mass; 1 mp = proton (hydrogen)
backgr_dens = 1.0e+15       -- background ion number density [cm^-3]
backgr_Z    = 1             -- background ion charge number (protons = 1)

-- Background bulk drift velocity expressed as a fraction of the Alfven speed.
-- Set all to 0 for a stationary (non-flowing) background.
-- Set backgr_flow_Ma_x = 8 to mimic a solar-wind-like super-Alfvenic flow.
backgr_flow_Ma_x = 0.0      -- drift along X in units of V_A0
backgr_flow_Ma_y = 0.0
backgr_flow_Ma_z = 0.0

-- ===========================================================================
-- EXTERNAL MAGNETIC FIELD
-- ===========================================================================
-- A uniform field along Z.  The Alfven speed V_A0 is derived from this.
Bz0 = 100.0                 -- background magnetic field [Gauss]

-- ===========================================================================
-- PLASMA CLOUD (expanding sphere)
-- ===========================================================================
-- The cloud represents an impulsively heated or laser-driven plasma blob.
-- Its ions expand radially with velocities from 0 to cloud_v_max = Ma * V_A0.
-- The Alfven Mach number Ma controls how strongly the cloud brakes.

cloud_mi       = 1 * cgs.mp -- cloud ion mass (same species as background here)
cloud_ions_num = 6.0e+18    -- total number of REAL ions in the cloud
cloud_Z        = 1          -- cloud ion charge number
Ma             = 2.0        -- Alfven Mach number (expansion speed / V_A0)

-- ===========================================================================
-- ALGORITHM SETTINGS
-- ===========================================================================
-- magnetic_field_alg: how we solve Faraday's law dB/dt = -curl(E)
--   FDTD -- finite-difference stencil, simple and local (default)
--   PSTD -- pseudo-spectral (FFT-based), no numerical dispersion but global
magnetic_field_alg = MagneticFieldAlg.FDTD

-- Verbosity level (controls console output):
--   0 = silent (no timing output at all)
--   1 = step markers only: "Step X of Y" on every step (default)
--   2 = verbose: step markers + intermediate timers (for debugging)
-- Set before calling require("sim_core").run():
--   pic_config.set_verbosity_level(0)  -- for quiet runs
--   pic_config.set_verbosity_level(2)  -- for detailed output

-- Isotropized curl (transversely-smoothed Cole-Karkkainen-style stencils in
-- the FDTD field solver; cancels the leading ANISOTROPIC term of the discrete
-- curl-curl dispersion error): OFF by default (C++ default reproduces the
-- classic 2-point Yee curl bit-for-bit). To enable in a case, call BEFORE
-- require("sim_core").run():
--   pic_config.set_isotropic_curl_enabled(true)

-- cold_electrons_enabled: if true, electron pressure and resistivity are zero.
-- Use this for initial testing or purely electromagnetic problems.
-- Set to false to enable Spitzer resistivity and electron heat conduction.
cold_electrons_enabled = true

-- Particle position jitter: adds a small random displacement to the regular
-- lattice to reduce artificial numerical noise from particle regularity.
-- Disabled here for reproducibility; enable for production runs.
cloud_jitter_enabled  = false
backgr_jitter_enabled = false

-- pstd_quick_test_steps: limit the run to N steps for a fast PSTD sanity check.
-- Set to nil for a full-length physical run.
pstd_quick_test_steps = nil

-- ===========================================================================
-- RESOLUTION PARAMETERS
-- ===========================================================================
-- These control the accuracy/cost trade-off.

-- points_on_RL: grid cells per ion Larmor radius.
-- Minimum ~5; use 10 for production accuracy.
points_on_RL = 10

-- h_R0: initial cloud radius expressed in grid steps.
-- Larger = better resolved cloud interior, but bigger grid.
h_R0 = 4.0

-- break_times: how many braking times T_b to simulate.
-- T_b = R_b / V_max where R_b is the magnetic braking radius.
-- 4 times is enough to see the cloud stop and form a cavity.
break_times = 4.0

-- Macro-particles per grid step (per axis).
-- The total number of macro-particles scales as parts_on_step^3 * volume.
-- More particles -> less statistical noise, but more memory and compute time.
cloud_parts_on_step  = 25   -- cloud:  25^3 * (4/3 pi h_R0^3) ~ 200k particles
backgr_parts_on_step = 4    -- background: 4^3 per cell * N_cells

-- ===========================================================================
-- LAUNCH SIMULATION
-- ===========================================================================
-- sim_core.run() calls all modules in the correct order:
--   physics.compute() -> init_particles.init() -> init_field.init() -> done
require("sim_core").run()

-- ===========================================================================
-- ELECTRON CLOSURE OVERRIDES (optional; must come AFTER run())
-- ===========================================================================
-- The electron-fluid closure has three tunable coefficients. Defaults are
-- assigned inside run() (physics.lua / sim_core.lua), so overrides only take
-- effect if set AFTER the run() call above. Full reference with measured
-- effects: sim/lib/README.md, section "Electron closure parameters";
-- textbook-style guide with derivations: docs/Eta_Chi_Kappa_Guide.md.
--
-- 1. eta -- resistivity [s]. Physical (Spitzer at Te_ref = 1 eV), computed
--    automatically; do not set `resistivity` by hand. It is the MAIN damper
--    of grid-scale field noise (diffusion number D = eta*c^2*tau/(4*pi*h^2)
--    ~ 0.04; the code caps eta at the explicit-stability limit D <= 0.1).
--    The Te-adaptive scaling eta ~ Te^-1.5 weakens damping at the hot front
--    and is the measured driver of the m=4 grid mode. Knobs:
--      pic_parameters.Spitzer_eta_floor_mult = 1.0
--        -- floor eta at base_eta; binds only where Te > Te_ref (the front).
--           Minimal-footprint fix for the m=4 grid mode (recommended).
--      pic_parameters.Spitzer_Te_ref = 0.0
--        -- constant eta everywhere (diagnostic switch). WARNING: the grid
--           dump normalizes its Te column by this value -> Te becomes inf.
--
-- 2. chi -- electron heat conduction [cm^2/s]:
--      pic_parameters.electron_thermal_conductivity = 1.93e+5
--    Physical Spitzer value for Te = 1 eV, n = 1e15 cm^-3. Implicit solver,
--    unconditionally stable at any value; 0 disables. NOTE: the effective
--    diffusivity is chi * clamp(Te/eta, 1e-6, 100) and at typical scales
--    the clamp saturates -> effectively 100*chi, uniform. Measured to have
--    NO effect on the m=4 grid mode under any eta model.
--
-- 3. kappa -- legacy-code-style interspecies friction (dimensionless):
--      pic_parameters.kappa_friction = 0.1   -- legacy M20D14_MATCH value
--      pic_parameters.kappa_B_ref    = Bz0   -- reference field B0 [G]
--    When both > 0 this REPLACES the Spitzer eta model with the legacy
--    closure: eta = kappa*B_ref/(c*e*ne) (no Te feedback, ~1/ne), friction
--    heating, and per-ion drag du/dt = kappa*Omega_ci*(Ue - u). The three
--    channels come as a package. Conversions: eta_equiv = kappa*B_ref/(c*e*n),
--    kappa_equiv = eta*c*e*n/B_ref (our Spitzer base_eta at n0 ~ kappa 16.5).
--    Keep Spitzer_Te_ref nonzero -- the Te dump column is normalized by it.
