-- main.lua -- WHISTLER_DISPERSION/FINE
-- =====================================================================
-- Purpose: docs/AHKASH_Comparison.md sec.5 -- does OpenPIC's leapfrog/
-- FDTD time scheme reproduce the parallel whistler dispersion relation
-- (docs/Hybrid_Condition_h_gt_cWpi.md, "Whistler-Branch Argument"):
--
--   omega = V_A * k * ( sqrt(1 + (k*d_i/2)^2) + k*d_i/2 )     ["+" branch]
--
-- FINE-regime companion of ../COARSE/main.lua -- identical physics
-- (n0, B0, k), identical box (length_x = 63.34197 cm, 7 wavelengths),
-- only h is halved (h = d_i instead of 2*d_i, the "resolving regime" of
-- docs/Hybrid_Condition_h_gt_cWpi.md, same regime as the HFINE
-- h_override cases in 09_DENSE_CLOUD). See ../COARSE/main.lua header for
-- the full derivation and analysis method (circular B seed selects the
-- whistler branch only; read Bx(t),By(t) from any one cell of the
-- existing all_z_<center>_grd_*.dat time series; the numerical rotation
-- period from Bx(t) zero-crossings is compared to the analytic T below
-- at BOTH resolutions -- if FINE does not track the analytic T better
-- than COARSE (or gets worse), that is evidence the time-integration
-- scheme itself struggles in the resolving regime, independent of the
-- m=4 front mode / ISOCURL findings in docs/M4_Remedies_Summary.md).
--
-- Derived parameters (same k, box as COARSE; only h differs):
--   d_i = 0.7200848 cm   V_A0 = 689757.08 cm/s   Wci = 957883.13 rad/s
--   k*d_i = 0.5 -> k = 0.6943627 rad/cm   lambda = 9.048852 cm
--   omega (whistler, "+" branch) = 613417.06 rad/s   T = 1.0242926e-05 s
--   h = d_i = 0.7197951 cm (rounded to fit 89 nodes exactly across the
--     same 63.34197 cm box as COARSE)
--   tau: Wci-limited (0.02/Wci) = 2.0879374e-08 s -- SAME as COARSE
--     (Wci does not depend on h) -> steps_per_period ~= 490.6
--   time_steps = 1500 (~3.06 periods), save_time_steps = 5
--
-- Run: opic main.lua   (from this directory)

package.path = package.path .. ";../../lib/?.lua;../../../../lib/?.lua"

local verify = require("verification_common")
local cgs = require("cgs")

local nodes = 89
local h = 0.7197950785863759
local B0 = 100.0
local n0 = 1.0e15
local mi = cgs.mp
local amp = 1.0e-6

local length_x = (nodes - 1) * h
local k = 2.0 * math.pi * 7.0 / length_x   -- same 7 wavelengths as COARSE

verify.run({
    case_name = "WHISTLER_DISPERSION_FINE: h = d_i, k*d_i = 0.5",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = n0,
    backgr_mi = mi,
    time_steps = 1500,
    save_time_steps = 5,
    save_whole_grid = true,

    field_fn = function(cell, x, y, z)
        cell.B.x = amp * B0 * math.cos(k * z)
        cell.B.y = amp * B0 * math.sin(k * z)
        cell.B.z = B0
    end,
})
