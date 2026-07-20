-- main.lua -- INFLOW_DENSITY_CONSERVATION
-- =====================================================================
-- Purpose: verify the X-axis particle-injection mechanism
-- (callbacks.lua: inject_background_particles/add_inflow_plane_x, wired
-- to on_iteration_begin) actually keeps the background density constant
-- once particles start leaving through the outflow wall. This is a pure
-- audit of an EXISTING mechanism, not a new one: the "backgr" group is
-- given a uniform drift velocity Vx (no external magnetic field, Bz0=0,
-- so there is no magnetic force at all -- the only thing moving the
-- particles is their initial drift), and the resulting flow is left
-- alone to run. Absent injection, the box would drain from the inflow
-- side and the interior density would fall monotonically. With
-- injection working correctly, a probe far from both walls should see
-- NP stay at backgr_dens throughout, because the plane of particles
-- leaving through the absorptive far wall is continuously replaced by
-- a plane entering through the near wall at the same rate (same density,
-- same velocity).
--
-- Setup: verification_common.lua builds a uniform "backgr"-only box (no
-- cloud, cloud_R0=0) with backgr_vx set directly (not derived from a
-- Mach number/Bz0, since Bz0=0 here) and automatically turns on
-- BoundaryKind.FlowBackground for "backgr" because backgr_v_abs > 0.
-- That, in turn, is what makes callbacks.lua's on_iteration_begin()
-- injection active (it only fires for backgr_vx ~= 0).
--
-- Derived numbers (h=1.0 cm, backgr_dens=1e15 cm^-3, backgr_vx=1e6 cm/s,
-- backgr_mi = mp, Bz0=0 so V_A0=0 and tau is set purely by the
-- advective CFL condition tau_cfl = 0.05*h/backgr_vx = 5e-8 s):
--   nodes = 21, length_x = 20 cm
--   box-crossing time = length_x/backgr_vx = 2e-5 s = 400 steps
--   time_steps = 2000 (5 box crossings -- enough to leave any initial
--     transient behind and see several steady-state injection cycles),
--   save_time_steps = 20 (100 snapshots)
--
-- Analysis: read the all_z_<center>_grd_*.dat equatorial slices (already
-- written every save step by diag.lua) and sample NP at the grid point
-- closest to the box center in X and Y (as far as possible from both
-- the inflow wall at x=0 and the absorptive outflow wall at x=length_x).
-- See tools/verify_inflow_density.py for the pass/fail criterion.
--
-- Run: opic main.lua   (from this directory)

package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")
local cgs = require("cgs")

verify.run({
    case_name = "INFLOW_DENSITY_CONSERVATION: X-axis injection keeps background density constant",
    nodes = 21,
    h = 1.0,
    Bz0 = 0.0,                 -- no external magnetic field, per request
    backgr_dens = 1.0e15,
    backgr_mi = cgs.mp,
    backgr_Z = 1,
    backgr_vx = 1.0e6,         -- uniform drift along x; only nonzero component
    U_scale = 1.0e6,           -- velocity scale: flow velocity (must be nonzero even when Bz0=0)
    E_scale = 1.0,             -- E scale: arbitrary nonzero (B=0 means no fields in this run)
    B_scale = 1.0,             -- B scale: arbitrary nonzero reference (actual B=0)
    time_steps = 2000,
    save_time_steps = 20,
    save_whole_grid = true,
    save_grid_z_plains = true,
})

-- verification_common.lua does not call opic_finalize_logging (it
-- bypasses sim_core.lua, so lua_params.txt is never written for
-- verification tasks) -- write the handful of values the check script
-- needs directly, from the globals verify.run() just set.
local f = io.open("test_params.txt", "w")
f:write("h = ", h, "\n")
f:write("length_x = ", length_x, "\n")
f:write("tau = ", pic_parameters.tau, "\n")
f:write("backgr_dens = ", backgr_dens, "\n")
f:write("backgr_vx = ", backgr_vx, "\n")
f:write("N_scale = ", pic_parameters.N_scale, "\n")
f:write("time_steps = ", pic_parameters.time_steps, "\n")
f:write("save_time_steps = ", pic_parameters.save_time_steps, "\n")
f:close()
