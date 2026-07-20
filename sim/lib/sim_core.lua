-- sim_core.lua -- Simulation initialization orchestrator
--
-- This module drives the startup sequence in the correct dependency order.
-- Think of it as the "main function" of the simulation setup.
--
-- ORDER MATTERS:
--   1. Register callbacks FIRST -- the engine may call them at any point.
--   2. Compute physics -- derives h, tau, grid size, particle counts.
--   3. Resize the grid -- required before any particle placement.
--   4. Check memory -- warn the user before a crash.
--   5. Place particles -- cloud then background.
--   6. Initialize fields -- B, E, Te, eta.
--   7. Set pic_parameters -- tell the C++ engine which algorithms to use.
--
-- Applicability checks and lua_params.txt logging are NOT done here --
-- see opic_finalize_logging() below: the C++ engine calls it automatically
-- after the WHOLE main.lua script (including anything the case does after
-- P.run() returns) has finished, so the log always reflects the final
-- pic_parameters actually used by the simulation.

local P = {}

function P.run(opts)
    -- opts: optional case-specific overrides (e.g. { h_override = 2.277 }),
    -- passed on to physics.compute(). No globals, no _G tables.
    opts = opts or {}
    require("print_mpi")

    -- proc_idx identifies this MPI rank (0 = master).
    -- print_mpi.print_root only prints from rank 0 to avoid duplicate output.
    proc_idx = pic_parameters.process_idx
    local pr = print_mpi.print_root
    local function line()
        pr(proc_idx, "----------------------------------------------------------------\n")
    end
    local function row(name, value, unit)
        pr(proc_idx, string.format("  %-22s %s%s\n", name, tostring(value), unit and (" " .. unit) or ""))
    end

    pr(proc_idx, "\nOpenPIC setup\n")
    line()

    -- Step 1: Register all engine callbacks before anything else.
    require("callbacks")

    -- Step 2: Derive all plasma physics, CFL, grid dimensions, particle counts.
    require("physics").compute(opts)

    -- Print the most important derived quantities for a quick sanity check.
    row("case", case_name or "(unnamed)")
    row("grid", string.format("%d x %d x %d", x_nodes, y_nodes, z_nodes))
    row("h", h, "cm")
    row("tau", tau, "s")
    row("time steps", pic_parameters.time_steps)
    row("B0", B0, "G")
    row("Ma", Ma)
    row("Delta", Delta)
    row("RL", RL, "cm")
    row("V_A0", V_A0, "cm/s")
    row("c/Wpi background", backgr_c_Wpi, "cm")
    row("c/Wpi cloud", cloud_c_Wpi, "cm")
    row("h/(c/Wpi bg)", h / backgr_c_Wpi)
    row("h/(c/Wpi cloud)", h / cloud_c_Wpi)
    row("CFL tau", CFL_tau, "s")
    row("Wci tau limit", Wci_tau, "s")
    row("local CFL tau", CFL_tau_loc, "s")
    row("cloud particles", cloud_parts_num)
    row("background particles", backgr_parts_num)

    -- Step 4: Memory estimate.
    -- Each particle: ~96 bytes.  Each grid cell: ~80 bytes.
    local mem_bytes = total_parts_num * 96 + x_nodes * y_nodes * z_nodes * 80
    row("estimated memory", math.ceil(mem_bytes/1024/1024), "MB")
    line()
    if mem_bytes > 8*1024*1024*1024 then
        pr(proc_idx, "WARNING: estimated memory may exceed 8 GB.\n")
        -- Only block for an answer in an interactive single-process run:
        -- under MPI no rank has a console to answer from, and a blocked
        -- io.read() would hang the whole job.
        if pic_parameters.process_num == 1 then
            pr(proc_idx, "Proceed? [y/n] ")
            local yn = io.read()
            if (yn == "n") or (yn == "N") then return end
        end
    end

    -- Step 5: Place particles (cloud sphere + background lattice).
    require("init_particles").init()

    -- Step 6: Initialize electromagnetic fields and electron temperature.
    pr(proc_idx, "\nInitializing fields\n")
    require("init_field").init()

    -- Step 7: Tell the C++ engine which numerical algorithms to use
    -- and set the physical normalization scales for output files.
    pic_parameters.tau              = tau
    pic_parameters.save_time_steps  = 10      -- save data every 10 steps
    pic_parameters.dens_cutoff      = 0.3 * backgr_dens  -- below this density -> skip physics
    pic_parameters.save_all_particles = false
    pic_parameters.save_whole_grid    = false
    pic_parameters.save_grid_x_plains = false -- save X=const planes?
    pic_parameters.save_grid_y_plains = false
    pic_parameters.save_grid_z_plains = true  -- save Z=const plane through cloud center

    -- Algorithm selection: push_method=Boris is the standard energy-conserving
    -- particle integrator. scatter_method selects one particle shape for BOTH
    -- deposition and field gather: Standard=CIC (default), NGP=nearest point,
    -- TSC=quadratic three-point B-spline.
    pic_parameters.push_method        = ParticlePushAlg.Boris
    pic_parameters.scatter_method     = opts.scatter_alg or ScatterAlg.Standard
    pic_parameters.magnetic_field_alg = magnetic_field_alg

    -- NOTE: Case-level CHOICES (isotropic_curl_enabled, Spitzer_eta_floor_mult,
    -- kappa_friction) are set by the case's main.lua via pic_config.set_*()
    -- calls -- the only sanctioned mechanism (no globals, no _G tables).
    -- Grid/physics overrides (h_override) arrive as run(opts).
    --
    -- resistivity/Spitzer_Te_ref are NOT case choices: physics.compute() just
    -- derived them from cold_electrons_enabled + Te_ref, and the case has no
    -- way to know them in advance. Only sim_core can forward them here.
    pic_parameters.resistivity    = resistivity
    pic_parameters.Spitzer_Te_ref = Spitzer_Te_ref

    -- Normalization scales: output files are written in these units.
    -- L_scale = c/Wpi (ion inertial length), T_scale = 1/Wci (cyclotron period).
    pic_parameters.L_scale = backgr_c_Wpi
    pic_parameters.T_scale = 1.0 / backgr_Wci
    pic_parameters.U_scale = V_A0
    pic_parameters.N_scale = backgr_dens
    pic_parameters.E_scale = (V_A0 / c) * B0
    pic_parameters.B_scale = B0

    pr(proc_idx, "\nSetup complete\n")
    line()
    io.flush()
end

-- Finalize logging: applicability checks + lua_params.txt. Called
-- AUTOMATICALLY by the C++ engine (bind_to_lua.cpp) once the ENTIRE
-- main.lua script has finished executing -- including any case-level
-- pic_parameters overrides that run AFTER P.run() returns (e.g.
-- Spitzer_eta_floor_mult, save_particle_diagnostics). This is the single
-- source of truth for "what did the simulation actually use": no case
-- file has to remember to re-log or re-check by hand, and no ordering
-- mistake inside sim_core.lua can silently log stale/default values
-- again (this exact class of bug -- resistivity/Spitzer_Te_ref showing
-- 0.0 in lua_params.txt despite the simulation using the real value --
-- happened twice this project before this was made automatic).
function _G.opic_finalize_logging()
    require("applicability").check()
    local f = io.open("lua_params.txt", "w")
    require("params_log").write(f)
    f:close()
end

return P
