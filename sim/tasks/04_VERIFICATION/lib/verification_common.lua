-- Shared helpers for small verification tasks.
--
-- These tasks intentionally bypass sim_core.lua: they build controlled
-- uniform-plasma setups directly so verification cases are not contaminated by
-- cloud initialization or production-case parameter derivation.

local P = {}

local cgs = require("cgs")
local physics = require("physics")

local function default(value, fallback)
    if value == nil then return fallback end
    return value
end

local function install_cylindrical_sources(params, nodes, h)
    local sources = params.cylindrical_sources or params.cylindrical_source
    if not sources then return end
    if sources.radius or sources.vx or sources.density then
        sources = { sources }
    end

    local base_on_iteration_begin = on_iteration_begin
    for _, source in ipairs(sources) do
        source.accum = 0.0
        source.group_name = source.group_name or "backgr"
        source.group_id   = pic_particle_groups:get_group(source.group_name).id
        source.part_dist = source.part_dist or h
        source.radius = source.radius or (4.0 * h)
        source.cy = default(source.cy, 0.5 * (nodes - 1) * h)
        source.cz = default(source.cz, 0.5 * (nodes - 1) * h)
        source.vx = default(source.vx, 0.0)
        source.vy = default(source.vy, 0.0)
        source.vz = default(source.vz, 0.0)
        source.density = default(source.density, params.backgr_dens)
        source.ni = source.density * h * h * h
        source.y0 = source.y0 or h
        source.y1 = source.y1 or ((nodes - 2) * h)
        source.z0 = source.z0 or h
        source.z1 = source.z1 or ((nodes - 2) * h)
    end

    function on_iteration_begin()
        if base_on_iteration_begin then base_on_iteration_begin() end
        if pic_parameters.process_idx ~= 0 then return end

        for _, source in ipairs(sources) do
            local vx_abs = math.abs(source.vx)
            if vx_abs > 0.0 and source.part_dist > 0.0 and source.ni > 0.0 then

                source.accum = source.accum + vx_abs * pic_parameters.tau
                if source.accum >= source.part_dist then

                    local planes = math.floor(source.accum / source.part_dist)
                    local residual = source.accum - planes * source.part_dist
                    local edge = source.edge_x or h
                    local x_in = (source.vx >= 0.0) and edge or ((nodes - 1) * h - edge)
                    local x_sign = (source.vx >= 0.0) and 1.0 or -1.0

                    for plane = (planes - 1), 0, -1 do
                        local x = x_in + x_sign * (residual + plane * source.part_dist)
                        for y = source.y0, source.y1, source.part_dist do
                        for z = source.z0, source.z1, source.part_dist do
                            local dy = y - source.cy
                            local dz = z - source.cz
                            if (dy * dy + dz * dz) <= source.radius * source.radius then
                                pic_particles:add(source.group_id, x, y, z,
                                    source.vx, source.vy, source.vz, source.ni)
                            end
                        end
                        end
                    end

                    source.accum = residual
                end
            end
        end
    end
end

function P.run(params)
    require("print_mpi")
    require("callbacks")

    proc_idx = pic_parameters.process_idx

    -- Particles below are created identically on every rank (no rank
    -- partitioning), so an MPI run would multiply density by process_num.
    if pic_parameters.process_num > 1 then
        error("verification tasks are single-process only: they fill pic_particles on every rank")
    end

    pi = math.pi
    c  = cgs.c
    e  = cgs.e
    me = cgs.me
    mp = cgs.mp

    case_name   = params.case_name or "verification"
    h           = default(params.h, 1.0)
    -- Per-axis node counts default to the single `nodes` value (cubic grid,
    -- the shape every existing verification case uses). A case may override
    -- one axis independently (e.g. to push an absorbing boundary further
    -- from a region of interest) via nodes_x/nodes_y/nodes_z without
    -- affecting the other two.
    local nodes   = default(params.nodes, 17)
    local nodes_x = default(params.nodes_x, nodes)
    local nodes_y = default(params.nodes_y, nodes)
    local nodes_z = default(params.nodes_z, nodes)
    local steps = default(params.time_steps, 4)
    local save_steps = default(params.save_time_steps, 1)

    backgr_mi   = default(params.backgr_mi, cgs.mp)
    backgr_dens = default(params.backgr_dens, 1.0)
    backgr_Z    = default(params.backgr_Z, 1.0)
    Bz0         = default(params.Bz0, 1.0)

    B0   = math.abs(Bz0)
    V_A0 = physics.alfven_speed(B0, backgr_dens, backgr_mi)
    backgr_Wci   = physics.ion_cyclotron_frequency(backgr_Z, B0, backgr_mi)
    backgr_Wpi   = physics.ion_plasma_frequency(backgr_Z, backgr_dens, backgr_mi)
    backgr_c_Wpi = c / backgr_Wpi

    backgr_vx = default(params.backgr_vx, 0.0)
    backgr_vy = default(params.backgr_vy, 0.0)
    backgr_vz = default(params.backgr_vz, 0.0)
    backgr_v  = DblVector(backgr_vx, backgr_vy, backgr_vz)
    backgr_v_abs = math.sqrt(backgr_vx * backgr_vx + backgr_vy * backgr_vy + backgr_vz * backgr_vz)

    backgr_Ex0 = -backgr_vy * Bz0 / c
    backgr_Ey0 =  backgr_vx * Bz0 / c
    backgr_Ez0 = 0.0

    pic_grid.step = h
    pic_grid:resize(nodes_x, nodes_y, nodes_z)
    pic_grid:set_boundary_state(CellState.cs_absorptive)

    length_x = pic_grid:length_x()
    length_y = pic_grid:length_y()
    length_z = pic_grid:length_z()

    cloud_x_node = math.floor(0.5 * nodes_x)
    cloud_y_node = math.floor(0.5 * nodes_y)
    cloud_z_node = math.floor(0.5 * nodes_z)
    cloud_x = h * cloud_x_node
    cloud_y = h * cloud_y_node
    cloud_z = h * cloud_z_node
    cloud_R0 = 0.0
    cloud_R0_Sq = 0.0
    cloud_W0 = 1.0

    backgr_flow_boundary_cells_x = 2
    backgr_flow_boundary_cells_yz = 1
    -- Injection plane phase (callbacks.lua: add_inflow_plane_x) must match the
    -- initial lattice phase below (x=(i+0.5)*h), or the two interleave with a
    -- 0.5h offset -- a density seam that then rigidly translates through the
    -- box at backgr_vx and shows up as a transient NP anomaly when it crosses
    -- a sampling point (diagnosed 2026-07-18 via INFLOW_DENSITY_CONSERVATION).
    -- Must stay >= h: check_particle.cpp's is_particle_can_scatter kills any
    -- particle whose home cell index is < kd=1, i.e. x < h -- 0.5*h (tried
    -- first) puts injected particles inside that kill zone and they vanish
    -- on the same step they're injected. 1.5*h is the smallest offset that
    -- is both phase-matched ((n+0.5)*h for n=1) and outside the kill zone.
    backgr_particle_boundary_width_x = 1.5 * h
    backgr_particle_boundary_width_y = h
    backgr_particle_boundary_width_z = h
    backgr_part_dist = h
    backgr_ni = backgr_dens * h * h * h

    local speed_limit = math.max(V_A0, backgr_v_abs, 1.0)
    local tau_cfl = 0.05 * h / speed_limit
    local tau_wci = 0.02 / math.max(backgr_Wci, 1.0e-30)
    local tau = default(params.tau, math.min(tau_cfl, tau_wci))
    pic_parameters.tau = tau
    pic_parameters.time_steps = steps
    pic_parameters.save_time_steps = save_steps
    pic_parameters.dens_cutoff = default(params.dens_cutoff, 0.1 * backgr_dens)
    pic_parameters.save_particle_diagnostics = default(params.save_particle_diagnostics, false)
    pic_parameters.save_all_particles = false
    pic_parameters.save_whole_grid = default(params.save_whole_grid, true)
    pic_parameters.save_grid_x_plains = default(params.save_grid_x_plains, false)
    pic_parameters.save_grid_y_plains = default(params.save_grid_y_plains, false)
    pic_parameters.save_grid_z_plains = default(params.save_grid_z_plains, true)
    pic_parameters.push_method = ParticlePushAlg.Boris
    pic_parameters.scatter_method = params.scatter_alg or ScatterAlg.Standard
    pic_parameters.magnetic_field_alg = params.magnetic_field_alg or MagneticFieldAlg.FDTD
    pic_parameters.resistivity = default(params.resistivity, 0.0)
    pic_parameters.Spitzer_Te_ref = default(params.Spitzer_Te_ref, 1.0)
    pic_parameters.electron_thermal_conductivity = 0.0

    pic_parameters.L_scale = default(params.L_scale, h)
    pic_parameters.T_scale = default(params.T_scale, 1.0 / math.max(backgr_Wci, 1.0e-30))
    pic_parameters.U_scale = default(params.U_scale, math.max(V_A0, 1.0))
    pic_parameters.N_scale = backgr_dens
    pic_parameters.E_scale = default(params.E_scale, math.max((V_A0 / c) * B0, 1.0))
    pic_parameters.B_scale = default(params.B_scale, math.max(B0, 1.0))

    pic_particle_groups:create_group("backgr", 1.0, 1.0, Diagnostics.save_grid_values)
    if backgr_v_abs > 0.0 then
        pic_particle_groups:set_boundary_kind("backgr", BoundaryKind.FlowBackground)
    end
    install_cylindrical_sources(params, nodes, h)

    -- Extra co-located ion populations (e.g. a drifting beam superposed on
    -- the "backgr" core), each filled with its own uniform density/velocity
    -- at the SAME lattice positions as the core -- params.extra_groups is a
    -- list of { name, density, velocity = DblVector(...), mass, charge,
    -- diagnostics }. Not spatially varying (unlike density_fn/velocity_fn
    -- on the core group) since every case needing this so far wants a
    -- uniform drifting population, not a spatial profile.
    local extra_groups = params.extra_groups or {}
    for _, eg in ipairs(extra_groups) do
        local charge = default(eg.charge, 1.0)
        local mass   = default(eg.mass, 1.0)
        pic_particle_groups:create_group(eg.name, charge, mass, eg.diagnostics or Diagnostics.save_grid_values)
        if eg.velocity and (math.abs(eg.velocity.x) + math.abs(eg.velocity.y) + math.abs(eg.velocity.z)) > 0.0 then
            pic_particle_groups:set_boundary_kind(eg.name, BoundaryKind.FlowBackground)
        end
    end

    local particle_positions = {}
    for i = 1, nodes_x - 2 do
    for j = 1, nodes_y - 2 do
    for k = 1, nodes_z - 2 do
        particle_positions[#particle_positions + 1] = {
            x = (i + 0.5) * h,
            y = (j + 0.5) * h,
            z = (k + 0.5) * h
        }
    end
    end
    end

    pic_particles.size = #particle_positions * (1 + #extra_groups)
    local backgr_group_id = pic_particle_groups:get_group("backgr").id
    for p = 1, #particle_positions do
        local pos = particle_positions[p]
        local v = backgr_v
        if params.velocity_fn then
            v = params.velocity_fn(pos.x, pos.y, pos.z)
        end
        local dens = backgr_dens
        if params.density_fn then
            dens = params.density_fn(pos.x, pos.y, pos.z)
        end
        pic_particles:set(p - 1, backgr_group_id, pos.x, pos.y, pos.z, v.x, v.y, v.z, dens * h * h * h)
    end

    local next_index = #particle_positions
    for gi, eg in ipairs(extra_groups) do
        local group_id = pic_particle_groups:get_group(eg.name).id
        for p = 1, #particle_positions do
            local pos = particle_positions[p]
            pic_particles:set(next_index, group_id, pos.x, pos.y, pos.z,
                eg.velocity.x, eg.velocity.y, eg.velocity.z, eg.density * h * h * h)
            next_index = next_index + 1
        end
    end

    for i = 0, nodes_x - 1 do
    for j = 0, nodes_y - 1 do
    for k = 0, nodes_z - 1 do
        local x = i * h
        local y = j * h
        local z = k * h
        local cell = pic_grid:at(i, j, k)
        local dens = backgr_dens
        if params.density_fn then
            dens = params.density_fn(x, y, z)
        end
        cell.NP = dens
        cell.UP = DblVector(backgr_vx, backgr_vy, backgr_vz)
        cell.UE = DblVector(backgr_vx, backgr_vy, backgr_vz)
        cell.E  = DblVector(backgr_Ex0, backgr_Ey0, backgr_Ez0)
        cell.B  = DblVector(0.0, 0.0, Bz0)
        cell.Te = pic_parameters.Spitzer_Te_ref
        cell.eta = pic_parameters.resistivity
        if params.field_fn then
            params.field_fn(cell, x, y, z)
        end
    end
    end
    end

    print_mpi.print_root(proc_idx, "\n*** verification task: ", case_name, " ***\n")
    print_mpi.print_root(proc_idx, "nodes=(", nodes_x, ",", nodes_y, ",", nodes_z, ") particles=", #particle_positions,
        " tau=", tau, " steps=", steps, "\n\n")
end

return P
