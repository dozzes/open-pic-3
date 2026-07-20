-- callbacks.lua -- Engine callback functions
--
-- The C++ engine calls these global Lua functions at specific points in the
-- simulation loop.  They are the main extension points for custom physics.
--
-- CALLBACK ORDER PER STEP
-- -----------------------
--   on_iteration_begin()          -- before field/particle update
--   ... (engine does B, E, push) ...
--   on_particles_moved_half_time()
--   ... (engine scatters density) ...
--   on_set_boundary_NP()          -- fix density at boundaries
--   on_set_boundary_UP()          -- fix ion velocity at boundaries
--   on_set_boundary_UE()          -- fix electron velocity
--   on_set_boundary_EF()          -- fix electric field
--   on_set_boundary_MF()          -- fix magnetic field
--   on_particles_moved_full_time()
--   on_iteration_end()            -- diagnostics, archiving

-- Enable callback tracing for debugging (set to true to see all callbacks)
local enable_callback_trace = false

-- Helper: print callback name if tracing is enabled
local function trace_callback(name)
    if enable_callback_trace then
        print_mpi.print_root(proc_idx, "\t\t" .. name .. "\n")
    end
end

-- ============================================================================
-- BACKGROUND INFLOW
-- ============================================================================
-- When the background has a bulk drift (backgr_flow_Ma_x != 0), particles
-- stream out of the right boundary and must be replenished from the left.
-- We accumulate how far the flow has moved since the last injection
-- and add new particle planes whenever the accumulation exceeds one spacing.

local inflow_accum   = 0.0   -- accumulated inflow displacement [cm]
local inflow_warned  = false  -- print the vy/vz warning only once
local inflow_injected_count = 0  -- track total planes injected (debug)

local function add_inflow_plane_x(x)
    -- Add one plane of background particles at position x (perpendicular to X).
    -- Group id and velocity components are resolved once outside the loop
    -- and passed to the raw-double add() overload (fast path, same as
    -- init_particles.lua -- see Analysis Notes 2026-07-11): avoids a
    -- group-name lookup and a DblVector construction per particle.
    local ey = backgr_particle_boundary_width_y or h
    local ez = backgr_particle_boundary_width_z or ey
    local backgr_group_id = pic_particle_groups:get_group("backgr").id
    -- Use backgr_vx/vy/vz directly; backgr_v may not be accessible
    local vx = backgr_vx or 0.0
    local vy = backgr_vy or 0.0
    local vz = backgr_vz or 0.0
    for y = ey, (length_y - ey), backgr_part_dist do
    for z = ez, (length_z - ez), backgr_part_dist do
        pic_particles:add(backgr_group_id, x, y, z, vx, vy, vz, backgr_ni)
    end
    end
end

local function inject_background_particles()
    if pic_parameters.process_idx ~= 0 then return end

    -- Guard: inflow only makes sense for an X-directed background flow.
    if not (backgr_part_dist and backgr_vx and backgr_ni) then return end
    if backgr_part_dist <= 0.0 or backgr_ni <= 0.0 then return end

    local vx_abs = math.abs(backgr_vx)
    if vx_abs <= 0.0 then return end

    -- Warn once if the flow has Y or Z components that we are not replenishing.
    if (backgr_vy ~= 0.0 or backgr_vz ~= 0.0) and not inflow_warned then
        print_mpi.print_root(proc_idx,
            "WARNING: inflow replenishment handles X-direction only.\n")
        inflow_warned = true
    end

    -- Accumulate displacement; inject whole planes when we have moved >= one spacing.
    inflow_accum = inflow_accum + vx_abs * pic_parameters.tau
    if inflow_accum < backgr_part_dist then return end

    local planes  = math.floor(inflow_accum / backgr_part_dist)
    local residual = inflow_accum - planes * backgr_part_dist

    -- Inject planes at the inflow boundary edge (+X for positive vx, -X for negative).
    local edge   = backgr_particle_boundary_width_x or h
    local x_in   = (backgr_vx >= 0.0) and edge or (length_x - edge)
    local x_sign = (backgr_vx >= 0.0) and 1.0 or -1.0

    -- A plane injected here (on_iteration_begin, before this step's push)
    -- still gets pushed by vx*tau later in this same step, same as every
    -- other particle. Without compensation it lands one step's displacement
    -- past the intended lattice phase, producing a 1-step (5% of one
    -- spacing here) phase offset from the pre-existing lattice; that small
    -- gap/pile-up then rigidly convects through the box and shows up as a
    -- transient NP deviation when it crosses a sampling point (diagnosed
    -- 2026-07-18 via INFLOW_DENSITY_CONSERVATION -- residual ~2.5% bump
    -- after the plane-vs-plane phase fix). Pre-shift backwards by one
    -- step's push so the post-push position matches the intended phase.
    local step_shift = x_sign * vx_abs * pic_parameters.tau

    for plane = (planes - 1), 0, -1 do
        add_inflow_plane_x(x_in + x_sign * (residual + plane * backgr_part_dist) - step_shift)
        inflow_injected_count = inflow_injected_count + 1
    end
    inflow_accum = residual
end

-- ============================================================================
-- ITERATION CALLBACKS
-- ============================================================================

function on_iteration_begin()
    -- Replenish background particles absorbed at the inflow boundary.
    inject_background_particles()
end

function on_iteration_end()
    local t  = pic_parameters.current_time_step
    local st = pic_parameters.save_time_steps

    -- Create the diagnostics directory on the very first step.
    if t == 1 then require("utils").ensure_dir("diag") end

    -- Run diagnostics on save steps and on the final step.
    if (t == 1) or (t % st == 0) or (t == st) then
        require("diag").on_step(t, st)
    end
end

function on_particles_moved_half_time()
    trace_callback("on_particles_moved_half_time")
end

function on_particles_moved_full_time()
    trace_callback("on_particles_moved_full_time")
end

-- ============================================================================
-- BOUNDARY CONDITION CALLBACKS
-- ============================================================================
-- These functions are called after each scatter/field-update pass.
-- They enforce the physical boundary conditions at the domain edges.
-- "nonperturbed" means: copy interior values into the boundary ghost cells
-- so that the boundary looks like an undisturbed plasma (open/absorbing BC).
-- "flow" variants also force the boundary cells to the background drift value.

function on_set_boundary_NP()
    trace_callback("on_set_boundary_NP")
    require("bound_cond").nonperturbed_flow_NP(pic_grid)
end

function on_set_boundary_group_NP(grid, group_name)
    trace_callback("on_set_boundary_group_NP")
    local bc    = require("bound_cond")
    local group = pic_particle_groups:get_group(group_name)
    if group.boundary_kind == BoundaryKind.FlowBackground then
        bc.nonperturbed_flow_NP(grid)
    else
        bc.nonperturbed_grp_NP(grid)
    end
end

function on_set_boundary_UP()
    trace_callback("on_set_boundary_UP")
    require("bound_cond").nonperturbed_flow_UP(pic_grid)
end

function on_set_boundary_group_UP(grid, group_name)
    trace_callback("on_set_boundary_group_UP")
    local bc    = require("bound_cond")
    local group = pic_particle_groups:get_group(group_name)
    if group.boundary_kind == BoundaryKind.FlowBackground then
        bc.nonperturbed_flow_UP(grid)
    else
        bc.nonperturbed_grp_UP(grid)
    end
end

function on_set_boundary_UE()
    trace_callback("on_set_boundary_UE")
    require("bound_cond").nonperturbed_UE(pic_grid)
end

function on_set_boundary_Te()
    trace_callback("on_set_boundary_Te")
    require("bound_cond").nonperturbed_Te(pic_grid)
end

function on_set_boundary_eta()
    trace_callback("on_set_boundary_eta")
    require("bound_cond").nonperturbed_eta(pic_grid)
end

function on_set_boundary_EF()
    trace_callback("on_set_boundary_EF")
    require("bound_cond").nonperturbed_EF(pic_grid)
end

function on_set_boundary_MF()
    trace_callback("on_set_boundary_MF")
    require("bound_cond").nonperturbed_MF(pic_grid)
end
