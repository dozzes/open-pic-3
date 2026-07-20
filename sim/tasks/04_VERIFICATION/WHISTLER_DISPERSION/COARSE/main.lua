-- main.lua -- WHISTLER_DISPERSION/COARSE
-- =====================================================================
-- Purpose: docs/AHKASH_Comparison.md sec.5 -- does OpenPIC's leapfrog/
-- FDTD time scheme reproduce the parallel whistler dispersion relation
-- (docs/Hybrid_Condition_h_gt_cWpi.md, "Whistler-Branch Argument"):
--
--   omega = V_A * k * ( sqrt(1 + (k*d_i/2)^2) + k*d_i/2 )     ["+" branch]
--
-- No cloud -- pure uniform background plasma + Bz0, seeded with a
-- single-k circularly polarized B perturbation (Bx,By) around Bz0. This
-- selects ONLY the right-hand/whistler branch (not a superposition with
-- the left-hand ion-cyclotron branch, which a linear Bx-only seed would
-- also excite) so a single probe point gives a clean rotation frequency,
-- not a beat pattern. Since the seed depends only on z, every cell in a
-- z=const plane has identical Bx,By at all times -- read any one cell
-- from the existing all_z_<center>_grd_*.dat time series (diag.lua) as
-- Bx(t), By(t) at the probe; extract the numerical rotation frequency
-- from successive zero-crossings of Bx(t) (period T_num = 2*pi/omega_num)
-- and compare to the analytic T above.
--
-- This is the COARSE-regime companion of ./../FINE/main.lua (h = d_i):
-- same physics (n0, B0, k), h = 2*d_i (the historical floor-pinned
-- default, physics.lua's h = max(2*c/Wpi, RL/points_on_RL)).
--
-- Derived parameters (n0=1e15 cm^-3, mi=mp, B0=100 G -- same convention
-- as 04_VERIFICATION/LINEAR_ALFVEN):
--   d_i = 0.7200848 cm   V_A0 = 689757.08 cm/s   Wci = 957883.13 rad/s
--   k*d_i = 0.5 -> k = 0.6943627 rad/cm   lambda = 9.048852 cm
--   omega (whistler, "+" branch) = 613417.06 rad/s   T = 1.0242926e-05 s
--   h = 2*d_i = 1.4395902 cm (rounded to fit 45 nodes exactly across
--     7 wavelengths: length_x = 63.34197 cm, nodes = 45)
--   tau: Wci-limited (0.02/Wci) = 2.0879374e-08 s -- same tau as FINE,
--     since Wci does not depend on h -> steps_per_period ~= 490.6
--   time_steps = 1500 (~3.06 periods), save_time_steps = 5
--     (~98 samples/period -- plenty for zero-crossing timing)
--
-- Run: opic main.lua   (from this directory)

package.path = package.path .. ";../../lib/?.lua;../../../../lib/?.lua"

local verify = require("verification_common")
local cgs = require("cgs")

local nodes = 45
local h = 1.4395901571727518
local B0 = 100.0
local n0 = 1.0e15
local mi = cgs.mp
local amp = 1.0e-6

local length_x = (nodes - 1) * h
local k = 2.0 * math.pi * 7.0 / length_x   -- 7 wavelengths across the box

verify.run({
    case_name = "WHISTLER_DISPERSION_COARSE: h = 2*d_i, k*d_i = 0.5",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = n0,
    backgr_mi = mi,
    time_steps = 1500,
    save_time_steps = 5,
    save_whole_grid = true,

    -- Circularly polarized seed around Bz0 -- selects the whistler
    -- ("+") branch only, no velocity_fn override (ions respond
    -- self-consistently through Ohm's law once the run starts).
    field_fn = function(cell, x, y, z)
        cell.B.x = amp * B0 * math.cos(k * z)
        cell.B.y = amp * B0 * math.sin(k * z)
        cell.B.z = B0
    end,
})
