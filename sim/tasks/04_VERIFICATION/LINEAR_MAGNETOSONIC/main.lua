package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")

local nodes = 33
local h = 1.0
local B0 = 100.0
local n0 = 1.0e15
local mi = require("cgs").mp
local amp = 1.0e-6
local Lx = (nodes - 1) * h
local k = 2.0 * math.pi / Lx
local vA = B0 / math.sqrt(4.0 * math.pi * n0 * mi)

verify.run({
    case_name = "LINEAR_MAGNETOSONIC",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = n0,
    backgr_mi = mi,
    time_steps = 8,
    save_time_steps = 1,
    save_whole_grid = true,

    velocity_fn = function(x, y, z)
        local dvx = amp * vA * math.sin(k * x)
        return DblVector(dvx, 0.0, 0.0)
    end,

    field_fn = function(cell, x, y, z)
        local phase = math.sin(k * x)
        cell.NP = n0 * (1.0 + amp * phase)
        cell.B.z = B0 * (1.0 + amp * phase)
    end,
})
