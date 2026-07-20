package.path = package.path .. ";../../lib/?.lua;../../../../lib/?.lua"

local verify = require("verification_common")
local cgs = require("cgs")

local nodes = 33
local h = 1.0
local B0 = 100.0
local n0 = 1.0e14
local mi = cgs.mp
local amp = 1.0e-6
local Lz = (nodes - 1) * h
local k = 2.0 * math.pi / Lz
local vA = B0 / math.sqrt(4.0 * math.pi * n0 * mi)
local cWpi = cgs.c / math.sqrt(4.0 * math.pi * n0 * cgs.e * cgs.e / mi)

print(string.format("Hybrid condition case: VIOLATED, h/(c/Wpi)=%.6g", h / cWpi))

verify.run({
    case_name = "HYBRID_CONDITION_VIOLATED_H_LT_CWPI",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = n0,
    backgr_mi = mi,
    time_steps = 40,
    save_time_steps = 10,
    save_whole_grid = true,

    velocity_fn = function(x, y, z)
        local dvx = amp * vA * math.sin(k * z)
        return DblVector(dvx, 0.0, 0.0)
    end,

    field_fn = function(cell, x, y, z)
        local dB = amp * B0 * math.sin(k * z)
        cell.B.x = dB
        cell.B.z = B0
    end,
})

