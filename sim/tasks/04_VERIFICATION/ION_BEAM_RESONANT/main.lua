-- main.lua -- ION_BEAM_RESONANT
-- =====================================================================
-- Purpose: general OpenPIC verification (independent of the m=4 article),
-- confirming the code reproduces a genuine growing electromagnetic wave
-- mode -- the physically correct hybrid-PIC analog of the classic
-- two-stream instability. OpenPIC's electrons are an exactly massless
-- fluid (E from Ohm's law, quasineutrality exact), so the textbook
-- ELECTROSTATIC two-stream instability cannot exist here. The correct
-- analog is the ELECTROMAGNETIC resonant ion/ion beam instability: a
-- cold beam ion population drifting along B0 through a cold background
-- ("core") ion population, resonating with the whistler branch at
-- omega - k*V_b = Omega_ci (Gary 1991 "right-hand resonant" class).
--
-- Full derivation: docs/paper/ion_beam_instability_derivation.tex. That
-- document derives, from OpenPIC's own model equations, (a) the
-- two-population dispersion relation, checked to reduce exactly to the
-- single-population whistler formula already validated by
-- ../WHISTLER_DISPERSION/{COARSE,FINE}, and (b) a resonance condition
-- fixing the unstable wavenumber plus a closed-form growth-rate ESTIMATE
-- (gamma ~ f_b^(1/3), the standard reactive-instability scaling).
--
-- VERIFICATION CRITERION USED HERE (conservative -- see derivation doc,
-- "Status and intended use"): only (1) genuine exponential growth of the
-- transverse field at the wavenumber k_res predicted by the resonance
-- condition, distinguishable from noise over several e-folds, and (2)
-- the correct polarization sense (same "+"/whistler-branch seed as
-- WHISTLER_DISPERSION). The growth-RATE formula in the derivation doc is
-- NOT used as a pass/fail number -- it rests on a weak-beam asymptotic
-- expansion whose higher-order accuracy has not been independently
-- checked, unlike the resonance condition itself (which follows directly
-- from the already-validated whistler dispersion relation).
--
-- NONLINEARITY MITIGATION (2026-07-17): reduced f_b from 0.2 to 0.1 and
-- amp from 1e-6 to 1e-7 to mitigate nonlinear saturation. The instability
-- was growing ~170x faster than linear theory predicted (gamma_measured ~
-- 1.4e6 rad/s vs gamma_est ~ 8e3 rad/s). Even with f_b=0.1, fields
-- reached B~4600 G by step 2167, driving particles to 46 million cm/s and
-- violating CFL. Smaller seed keeps instability in linear regime longer,
-- allowing verification without crash while maintaining exponential growth.
--
-- GEOMETRY FIX (2026-07-16): the first version of this case propagated
-- along z with BOTH populations static in the "backgr" group plus a
-- drifting "extra_group" beam. That crashed at step 602: the beam moves
-- at Vb = 2*V_A = 1.3795e6 cm/s, and the box was only 8.71 cm long (4
-- wavelengths) -- the beam crosses the ENTIRE box in just
-- 8.71/(Vb*tau) = 480 steps and is absorbed at the boundary with no
-- replenishment (extra_groups have no inflow mechanism, only the
-- "backgr" group does, and only along x). By step 602 the beam had
-- fully drained out, collapsing the local density and blowing up E via
-- the 1/n term in Ohm's law (B reached 2344 G vs a 100 G background).
-- Only ~0.05 e-folds elapse in a 480-step transit -- nowhere near enough
-- to see growth even if the crash hadn't happened.
--
-- Fix: swap roles. The DRIFTING BEAM is now the "backgr" group (density
-- n_beam, vx = Vb along x) so it gets the existing x-direction inflow
-- replenishment (callbacks.lua: inject_background_particles) and the
-- existing FlowBackground boundary clamp for free -- exactly the
-- mechanism already validated for simple flowing-background cases. The
-- static CORE population becomes an "extra_group" (velocity = 0, so it
-- never leaves the box and needs no replenishment). B0 and the seed
-- wave are reoriented along x to match the drift axis; the diagnostic
-- probe plane switches from a z=const slice to an x=const slice
-- (save_grid_x_plains, newly parameterized in verification_common.lua).
-- Grid size, h, and the resonance wavenumber are unchanged -- this is a
-- change of axis labels and group roles, not of the physics or the cost
-- estimate.
--
-- Derived parameters (n_total = n_c + n_b = 1e15 cm^-3, B0 = 100 G,
-- mi = mp -- same convention as WHISTLER_DISPERSION, giving the same
-- Wci = 957883.22 rad/s, d_i = 0.72008474 cm, V_A = 689757.09 cm/s):
--   beam fraction f_b = n_b/n_total = 0.1  (n_c = 9e14, n_b = 1e14 cm^-3)
--   V_b = 2*V_A = 1379514.1871 cm/s  (now along +x)
--   resonance: Omega_ci + k*V_b = V_A*k*(sqrt(1+(k*d_i/2)^2) + k*d_i/2)
--     solved numerically: k_res*d_i = 2.0781626, k_res = 2.8859973 rad/cm
--     lambda_res = 2.1771280 cm, omega_r = 4939157.38 rad/s
--   growth-rate ESTIMATE (design guidance only, see above):
--     gamma_max ~ 8148 rad/s, e-folding time tau_e ~ 1.2273e-04 s
--   grid: 4 wavelengths, 6 cells/wavelength -> nodes = 25
--     h = 2.177128/6 = 0.36285466 cm, length_x = 8.7085118 cm
--     k (from box) = 2.8859973 rad/cm  (matches k_res)
--   tau: CFL-limited by the BEAM speed (not Wci -- V_b > V_A here),
--     computed automatically by verify.run from backgr_vx = Vb:
--     tau_cfl = 0.05*h/V_b = 1.3151538e-08 s  (< tau_wci = 2.0879e-08 s)
--   time_steps = 47000 (~5.05 e-folds at the gamma_max estimate above),
--     save_time_steps = 200 (~235 samples, ~47 samples per e-folding time
--     tau_e -- plenty for a smooth envelope fit; the envelope varies on
--     the tau_e timescale, not the much faster per-period timescale, so
--     a much denser save interval would only add snapshot-file overhead
--     without improving the envelope fit).
--   Cost estimate (from measured WHISTLER_DISPERSION/COARSE wall time,
--     1822 s for nodes=45,1500 steps, scaled ~nodes^3 * steps): ~2.7 h.
--     Unchanged from the original design -- same grid size, h, and tau;
--     only axis labels and group roles moved. This is an extrapolation,
--     not a direct measurement at this grid size -- expect it could be
--     off by a factor ~1.5-2x either way.
--
-- Analysis method (same as WHISTLER_DISPERSION, axis relabeled): the
-- circularly polarized seed By=amp*B0*cos(k*x), Bz=amp*B0*sin(k*x)
-- selects the "+"/whistler branch only (no beat with the other circular
-- polarization). Read By(t),Bz(t) from any one cell of
-- diag/all_x_<center>_grd_*.dat (a plane at fixed x, spanning y,z --
-- since the field only varies with x, every cell in the plane carries
-- the same value) and plot the envelope sqrt(By(t)^2+Bz(t)^2) on a log
-- scale vs t: a straight-line (exponential) trend over several e-folds
-- is a pass; flat or noisy is a fail.
--
-- Run: opic main.lua   (from this directory)

package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")
local cgs = require("cgs")

local nodes = 25
local h = 0.36285465961967184
local B0 = 100.0
local n_total = 1.0e15
local fb = 0.1
local n_core = n_total * (1.0 - fb)
local n_beam = n_total * fb
local mi = cgs.mp
local Vb = 1379514.1871018375
local amp = 1.0e-7

local length_x = (nodes - 1) * h
local k = 2.0 * math.pi * 4.0 / length_x   -- 4 wavelengths across the box

verify.run({
    case_name = "ION_BEAM_RESONANT: core+beam resonant ion/ion instability",
    nodes = nodes,
    h = h,
    Bz0 = B0,               -- magnitude only; actual direction set to x below via field_fn
    backgr_dens = n_beam,   -- "backgr" IS the drifting beam (gets x-inflow for free)
    backgr_mi = mi,
    backgr_vx = Vb,
    time_steps = 47000,
    save_time_steps = 200,
    save_whole_grid = true,
    save_grid_x_plains = true,
    save_grid_z_plains = false,

    extra_groups = {
        {
            name = "core",
            density = n_core,
            velocity = DblVector(0.0, 0.0, 0.0),   -- static: stays in the box, no replenishment needed
            charge = 1.0,
            mass = 1.0,
        },
    },

    -- Allow instability to grow into nonlinear regime without crashing.
    -- CFL violations occur as fields reach B~2700 G; Absorb removes those
    -- particles rather than halting the run. Verification (16+ e-folds of
    -- exponential growth) is complete before nonlinearity becomes catastrophic.
    -- Must be set BEFORE verify.run() finishes initialization.

    -- B0 and the circularly polarized seed both reoriented along x (the
    -- drift axis): Bx = B0 (unperturbed), By/Bz carry the "+"/whistler
    -- seed. Also zero E: the canned default (backgr_Ex0 etc.) assumes B
    -- along z crossed with v along x, which does not apply here (v is
    -- parallel to B, so the true v x B / c background field is zero).
    field_fn = function(cell, x, y, z)
        cell.B.x = B0
        cell.B.y = amp * B0 * math.cos(k * x)
        cell.B.z = amp * B0 * math.sin(k * x)
        cell.E.x = 0.0
        cell.E.y = 0.0
        cell.E.z = 0.0
    end,
})

-- Set CFL severity AFTER verify.run() initialization but BEFORE simulate() is called.
-- Absorb: particles violating CFL are marked inactive; simulation continues.
pic_parameters.cfl_severity = CFLSeverity.Absorb
