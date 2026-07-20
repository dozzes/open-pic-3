package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")

local nodes = 33
local h = 1.0
local B0 = 100.0
local amp = 1.0e-6
local L = (nodes - 1) * h
local k = 2.0 * math.pi / L

verify.run({
    case_name = "ORSZAG_TANG_2D",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = 1.0e15,
    time_steps = 8,
    save_time_steps = 1,
    save_whole_grid = true,

    velocity_fn = function(x, y, z)
        return DblVector(-amp * math.sin(k * y), amp * math.sin(k * x), 0.0)
    end,

    field_fn = function(cell, x, y, z)
        cell.B.x = -amp * B0 * math.sin(k * y)
        cell.B.y =  amp * B0 * math.sin(2.0 * k * x)
        cell.B.z = B0
    end,
})
