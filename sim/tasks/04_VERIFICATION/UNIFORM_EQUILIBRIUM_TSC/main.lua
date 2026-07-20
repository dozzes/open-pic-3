package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")

-- Runtime consistency check for the complete TSC gather/deposit path.
-- A spatially uniform equilibrium must remain uniform for any normalized
-- particle shape.
verify.run({
    case_name = "UNIFORM_EQUILIBRIUM_TSC",
    nodes = 17,
    h = 1.0,
    Bz0 = 1.0,
    backgr_dens = 1.0,
    time_steps = 4,
    save_time_steps = 1,
    save_whole_grid = true,
    scatter_alg = ScatterAlg.TSC,
})
