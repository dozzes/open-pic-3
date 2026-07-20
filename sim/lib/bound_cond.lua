--[[
NAME: bound_cond
NOTES: Boundary conditions optimized for Staggered C-Grid topology.
]]

local P = {}

local function for_boundary_layers(grid, fn)
    local sx, sy, sz = grid:size_x() - 1, grid:size_y() - 1, grid:size_z() - 1
    local x_width = rawget(_G, "backgr_flow_boundary_cells_x") or 2
    local yz_width = rawget(_G, "backgr_flow_boundary_cells_yz") or 1

    for i = 0, sx do
    for j = 0, sy do
    for k = 0, sz do
        if i < x_width or i > sx - x_width or
           j < yz_width or j > sy - yz_width or
           k < yz_width or k > sz - yz_width then
            fn(grid:at(i, j, k))
        end
    end
    end
    end
end

-- 1. Helper for Cell-Centered values (NP, Te, eta)
-- Live at (i + 1/2, j + 1/2, k + 1/2)
--
-- Each of the three passes (Z, Y, X) clamps ALL THREE axis indices to the
-- true interior, not just its own axis. At a face cell this reduces to a
-- plain same-axis extension (the other two clamps are no-ops), but at an
-- edge or corner the source is the fully-diagonal interior neighbor -- the
-- physically correct extension of a d/dn=0 condition along every outward
-- normal at once. Since every clamped source index always lands strictly
-- inside [2,sx-3]x[2,sy-3]x[2,sz-3] regardless of which axis's ghost zone
-- (i,j,k) is in, no pass ever reads a cell another pass has written -- so
-- there is no cross-pass read-after-write hazard, and all three passes agree
-- on the same value for cells they both touch (found 2026-07-20 diagnosing a
-- domain-edge field blow-up in a finite-Te run: same-axis-only extension let
-- edge/corner cells inherit an already-boundary-forced neighbor instead of
-- genuine interior physics).
local function sync_centered(grid, field)
    local sx, sy, sz = grid:size_x() - 1, grid:size_y() - 1, grid:size_z() - 1

    local function clampX(idx)
        if idx < 2 then return 2 end
        if idx > sx - 3 then return sx - 3 end
        return idx
    end
    local function clampY(idx)
        if idx < 2 then return 2 end
        if idx > sy - 3 then return sy - 3 end
        return idx
    end
    local function clampZ(idx)
        if idx < 2 then return 2 end
        if idx > sz - 3 then return sz - 3 end
        return idx
    end

    -- Z
    for i = 0, sx do
        local si = clampX(i)
        for j = 0, sy do
            local sj = clampY(j)
            local v_top = grid:at(si, sj, sz-3)[field]
            grid:at(i,j,sz-2)[field], grid:at(i,j,sz-1)[field], grid:at(i,j,sz)[field] = v_top, v_top, v_top
            local v_bot = grid:at(si, sj, 2)[field]
            grid:at(i,j,0)[field], grid:at(i,j,1)[field] = v_bot, v_bot
        end
    end
    -- Y
    for i = 0, sx do
        local si = clampX(i)
        for k = 0, sz do
            local sk = clampZ(k)
            local v_top = grid:at(si, sy-3, sk)[field]
            grid:at(i,sy-2,k)[field], grid:at(i,sy-1,k)[field], grid:at(i,sy,k)[field] = v_top, v_top, v_top
            local v_bot = grid:at(si, 2, sk)[field]
            grid:at(i,1,k)[field], grid:at(i,0,k)[field] = v_bot, v_bot
        end
    end
    -- X
    for j = 0, sy do
        local sj = clampY(j)
        for k = 0, sz do
            local sk = clampZ(k)
            local v_top = grid:at(sx-3, sj, sk)[field]
            grid:at(sx-2, j, k)[field], grid:at(sx-1, j, k)[field], grid:at(sx, j, k)[field] = v_top, v_top, v_top
            local v_bot = grid:at(2, sj, sk)[field]
            grid:at(1, j, k)[field], grid:at(0, j, k)[field] = v_bot, v_bot
        end
    end
end

-- 2. Helper for Face/Edge components (Ex, UPx, Bx etc.)
-- This handles the offset correctly: if we are on Face X, the X-boundary needs special care
--
-- Same diagonal-clamp design as sync_centered above: each pass clamps all
-- three axis indices to the interior (clampX/Y/Z), not just its own axis, so
-- edge/corner cells draw from the true diagonal interior neighbor. The write
-- zone here is symmetric (2 cells each side, indices {0,1} and {s-1,s}); only
-- the source DEPTH (offX/offY/offZ, 2 or 3 depending on which component is
-- being synced) differs from sync_centered. Source indices always land in
-- [offX,sx-offX] etc., strictly outside {0,1,s-1,s} on every axis, so -- as
-- above -- no pass ever reads a cell any pass has written, regardless of
-- Z/Y/X order.
local function sync_component(grid, field_name, comp)
    local sx, sy, sz = grid:size_x() - 1, grid:size_y() - 1, grid:size_z() - 1

    -- Logic: if we sync 'x' component, the X-loop offset is -2 (since it's on the face),
    -- while Y and Z loops remain -3.
    local offX = (comp == "x") and 2 or 3
    local offY = (comp == "y") and 2 or 3
    local offZ = (comp == "z") and 2 or 3

    local function clampX(idx)
        if idx < 2 then return offX end
        if idx > sx - 2 then return sx - offX end
        return idx
    end
    local function clampY(idx)
        if idx < 2 then return offY end
        if idx > sy - 2 then return sy - offY end
        return idx
    end
    local function clampZ(idx)
        if idx < 2 then return offZ end
        if idx > sz - 2 then return sz - offZ end
        return idx
    end

    -- Z-boundary
    for i = 0, sx do
        local si = clampX(i)
        for j = 0, sy do
            local sj = clampY(j)
            local v = grid:at(si, sj, sz - offZ)[field_name][comp]
            grid:at(i,j,sz-1)[field_name][comp], grid:at(i,j,sz)[field_name][comp] = v, v
            v = grid:at(si, sj, offZ)[field_name][comp]
            grid:at(i,j,0)[field_name][comp], grid:at(i,j,1)[field_name][comp] = v, v
        end
    end

    -- Y-boundary
    for i = 0, sx do
        local si = clampX(i)
        for k = 0, sz do
            local sk = clampZ(k)
            local v = grid:at(si, sy - offY, sk)[field_name][comp]
            grid:at(i,sy-1,k)[field_name][comp], grid:at(i,sy,k)[field_name][comp] = v, v
            v = grid:at(si, offY, sk)[field_name][comp]
            grid:at(i,0,k)[field_name][comp], grid:at(i,1,k)[field_name][comp] = v, v
        end
    end

    -- X-boundary
    for j = 0, sy do
        local sj = clampY(j)
        for k = 0, sz do
            local sk = clampZ(k)
            local v = grid:at(sx - offX, sj, sk)[field_name][comp]
            grid:at(sx-1, j, k)[field_name][comp], grid:at(sx, j, k)[field_name][comp] = v, v
            v = grid:at(offX, sj, sk)[field_name][comp]
            grid:at(0, j, k)[field_name][comp], grid:at(1, j, k)[field_name][comp] = v, v
        end
    end
end

-- Public Interface
function P.nonperturbed_NP(grid) sync_centered(grid, "NP") end
function P.nonperturbed_Te(grid) sync_centered(grid, "Te") end
function P.nonperturbed_eta(grid) sync_centered(grid, "eta") end

function P.nonperturbed_flow_NP(grid)
    P.nonperturbed_NP(grid)

    local dens = rawget(_G, "backgr_dens")
    if dens == nil then return end

    for_boundary_layers(grid, function(cell)
        cell.NP = dens
    end)
end

function P.nonperturbed_UP(grid)
    sync_component(grid, "UP", "x")
    sync_component(grid, "UP", "y")
    sync_component(grid, "UP", "z")
end

function P.nonperturbed_flow_UP(grid)
    P.nonperturbed_UP(grid)

    local vx = rawget(_G, "backgr_vx") or 0.0
    local vy = rawget(_G, "backgr_vy") or 0.0
    local vz = rawget(_G, "backgr_vz") or 0.0
    if vx == 0.0 and vy == 0.0 and vz == 0.0 then return end

    for_boundary_layers(grid, function(cell)
        cell.UP.x = vx
        cell.UP.y = vy
        cell.UP.z = vz
    end)
end

function P.nonperturbed_UE(grid)
    sync_component(grid, "UE", "x")
    sync_component(grid, "UE", "y")
    sync_component(grid, "UE", "z")
end

function P.nonperturbed_MF(grid)
    sync_component(grid, "B", "x")
    sync_component(grid, "B", "y")
    sync_component(grid, "B", "z")
end

function P.nonperturbed_EF(grid)
    sync_component(grid, "E", "x")
    sync_component(grid, "E", "y")
    sync_component(grid, "E", "z")
end

-- Aliases for groups
P.nonperturbed_grp_NP = P.nonperturbed_NP
P.nonperturbed_grp_UP = P.nonperturbed_UP

return P
