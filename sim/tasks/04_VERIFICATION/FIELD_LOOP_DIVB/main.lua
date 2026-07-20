package.path = package.path .. ";../lib/?.lua;../../../lib/?.lua"

local verify = require("verification_common")

local nodes = 33
local h = 1.0
local B0 = 100.0
local amp = 1.0e-6
local cx = 0.5 * (nodes - 1) * h
local cy = cx
local radius = 0.22 * (nodes - 1) * h

verify.run({
    case_name = "FIELD_LOOP_DIVB",
    nodes = nodes,
    h = h,
    Bz0 = B0,
    backgr_dens = 1.0e15,
    backgr_vx = 1.0e3,
    backgr_vy = 5.0e2,
    time_steps = 8,
    save_time_steps = 1,
    save_whole_grid = true,

    field_fn = function(cell, x, y, z)
        local dx = x - cx
        local dy = y - cy
        local r = math.sqrt(dx * dx + dy * dy)
        if r > 1.0e-12 and r < radius then
            local a = amp * (radius - r) / radius
            cell.B.x = cell.B.x - a * dy / r
            cell.B.y = cell.B.y + a * dx / r
        end
    end,
})
