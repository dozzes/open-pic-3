package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")

verify.run({
    case_name = "UNIFORM_EQUILIBRIUM",
    nodes = 17,
    h = 1.0,
    Bz0 = 1.0,
    backgr_dens = 1.0,
    time_steps = 4,
    save_time_steps = 1,
    save_whole_grid = true,
})
