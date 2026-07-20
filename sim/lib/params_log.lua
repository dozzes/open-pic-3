-- params_log.lua -- Parameter logging to lua_params.txt
--
-- WHY LOG PARAMETERS?
-- -------------------
-- When you look at a diag/ folder six months later you need to know exactly
-- what h, tau, particle count, and algorithm were used.  This file writes
-- every derived quantity to lua_params.txt in the run directory.
--
-- Call write(f) AFTER physics.compute() and particle counting are done.

local P = {}

function P.write(f)
    local w = function(...) f:write(...) end
    local backgr_Wpe = e * 2 * math.sqrt(pi * backgr_dens / me)

    -- Case identification
    w("case_name                        = ", case_name, "\n")
    w("magnetic_field_alg               = ",
      (magnetic_field_alg == MagneticFieldAlg.PSTD) and "PSTD" or "FDTD", "\n")
    w("particle_shape                   = ",
      (pic_parameters.scatter_method == ScatterAlg.NGP) and "NGP"
      or ((pic_parameters.scatter_method == ScatterAlg.TSC) and "TSC" or "CIC"), "\n")
    w("\n")

    -- Lua configuration (from main.lua globals before run())
    w("--- Lua Configuration (from main.lua) ---\n")
    w("isotropic_curl_enabled           = ", tostring(isotropic_curl_enabled or false), "\n")
    w("cold_electrons_enabled           = ", tostring(cold_electrons_enabled), "\n")
    w("cloud_jitter_enabled             = ", tostring(cloud_jitter_enabled), "\n")
    w("backgr_jitter_enabled            = ", tostring(backgr_jitter_enabled), "\n")
    w("pstd_quick_test_steps            = ", tostring(pstd_quick_test_steps), "\n")
    w("\n")

    -- Physics/algorithm parameters (set via pic_config.set_*() before run())
    w("--- Parameters from pic_config (applied by C++) ---\n")
    w("resistivity                      = ", tostring(pic_parameters.resistivity), "\n")
    w("Spitzer_Te_ref                   = ", tostring(pic_parameters.Spitzer_Te_ref), "  K\n")
    w("Spitzer_eta_floor_mult           = ", tostring(pic_parameters.Spitzer_eta_floor_mult), "\n")
    w("electron_thermal_conductivity    = ", tostring(pic_parameters.electron_thermal_conductivity), "\n")
    w("hyper_resistivity                = ", tostring(pic_parameters.hyper_resistivity), "\n")
    w("kappa_friction                   = ", tostring(pic_parameters.kappa_friction), "\n")
    w("kappa_B_ref                      = ", tostring(pic_parameters.kappa_B_ref), "\n")
    w("random_seed                      = ", tostring(pic_parameters.random_seed), "\n")
    w("\n")

    -- Grid and time step
    w("h                = ", h,
      "  (must be > ", cloud_c_Wpi, ", backgr_c_Wpi = ", backgr_c_Wpi, ")\n")
    w("tau              = ", tau, "\n")
    w("CFL_tau          = ", CFL_tau,     "  (FDTD wave-speed limit)\n")
    w("PSTD_tau         = ", PSTD_tau,    "  (Alfven CFL for PSTD)\n")
    w("Wci_tau          = ", Wci_tau,     "  (cyclotron period / 10)\n")
    w("CFL_tau_loc      = ", CFL_tau_loc, "  (dispersive whistler limit)\n")
    w("\n")

    -- Grid dimensions
    w("box_break_times = ", box_break_times, "  (break_times used for box sizing only)\n")
    w("x_nodes = ", x_nodes, "\n")
    w("y_nodes = ", y_nodes, "\n")
    w("z_nodes = ", z_nodes, "\n")
    w("length_x = ", length_x, "  cm\n")
    w("length_y = ", length_y, "  cm\n")
    w("length_z = ", length_z, "  cm\n")
    w("cloud_x = ", cloud_x, "  cm\n")
    w("cloud_y = ", cloud_y, "  cm\n")
    w("cloud_z = ", cloud_z, "  cm\n")
    w("cloud_center_shift_x_h = ", cloud_center_shift_x_h, "\n")
    w("cloud_center_shift_y_h = ", cloud_center_shift_y_h, "\n")
    w("cloud_center_shift_z_h = ", cloud_center_shift_z_h, "\n")
    w("cloud_x_node = ", cloud_x_node, "\n")
    w("cloud_y_node = ", cloud_y_node, "\n")
    w("cloud_z_node = ", cloud_z_node, "\n")
    w("\n")

    -- Background plasma
    w("backgr_dens      = ", backgr_dens,    "  cm^-3\n")
    w("backgr_Wpi       = ", backgr_Wpi,     "  rad/s\n")
    w("backgr_c_Wpi     = ", backgr_c_Wpi,   "  cm  (ion inertial length)\n")
    w("backgr_Wci       = ", backgr_Wci,     "  rad/s  (cyclotron frequency)\n")
    w("backgr_Wci*tau   = ", backgr_Wci*tau, "  (must be < 0.1)\n")
    w("backgr_Wpe       = ", backgr_Wpe,     "  rad/s  (electron plasma freq)\n")
    w("backgr_c_Wpe     = ", c/backgr_Wpe,   "  cm  (electron skin depth)\n")
    w("backgr_ni        = ", backgr_ni, "\n")
    w("backgr_v_abs     = ", backgr_v_abs, "  cm/s\n")
    w("backgr_vx        = ", backgr_vx, "\n")
    w("backgr_vy        = ", backgr_vy, "\n")
    w("backgr_vz        = ", backgr_vz, "\n")
    w("backgr_Ex0       = ", backgr_Ex0, "\n")
    w("backgr_Ey0       = ", backgr_Ey0, "\n")
    w("backgr_Ez0       = ", backgr_Ez0, "\n")
    w("backgr_inflow_x  = ", tostring(backgr_vx ~= 0.0), "\n")
    w("backgr_flow_Ma_x = ", backgr_flow_Ma_x, "\n")
    w("backgr_flow_Ma_y = ", backgr_flow_Ma_y, "\n")
    w("backgr_flow_Ma_z = ", backgr_flow_Ma_z, "\n")
    w("backgr_flow_boundary_cells_x     = ", backgr_flow_boundary_cells_x, "\n")
    w("backgr_flow_boundary_cells_yz    = ", backgr_flow_boundary_cells_yz, "\n")
    w("backgr_particle_boundary_cells_x = ", backgr_particle_boundary_cells_x, "\n")
    w("backgr_particle_boundary_cells_yz= ", backgr_particle_boundary_cells_yz, "\n")
    w("\n")

    -- Cloud plasma
    w("cloud_dens       = ", cloud_dens,    "  cm^-3\n")
    w("cloud_Wpi        = ", cloud_Wpi,     "  rad/s\n")
    w("cloud_c_Wpi      = ", cloud_c_Wpi,   "  cm\n")
    w("cloud_Wci        = ", cloud_Wci,     "  rad/s\n")
    w("cloud_Wci*tau    = ", cloud_Wci*tau, "  (must be < 0.2)\n")
    w("cloud_c_Wci      = ", c/cloud_Wci,   "  cm\n")
    w("cloud_R0         = ", cloud_R0,      "  cm\n")
    w("cloud_W0         = ", cloud_W0,      "  erg  (reference kinetic energy)\n")
    w("cloud_part_dist  = ", cloud_part_dist, "\n")
    w("cloud_parts_num  = ", cloud_parts_num, "\n")
    w("cloud_ni         = ", cloud_ni,        "\n")
    w("\n")

    -- Characteristic radii and times
    w("RL     = ", RL,  "  cm  (cloud Larmor radius)\n")
    w("RL_b   = ", RL_b,"  cm  (background Larmor radius)\n")
    w("Rb     = ", Rb,  "  cm  (magnetic braking radius)\n")
    w("Rb1    = ", Rb1, "  cm  (mass-loading braking radius)\n")
    w("Rg     = ", Rg,  "  cm  (gas-dynamic braking radius)\n")
    w("Rt     = ", Rt,  "  cm  (actual braking radius = min(Rb,Rg))\n")
    w("Tb     = ", Tb,  "  s   (magnetic braking time)\n")
    w("Tg     = ", Tg,  "  s   (gas-dynamic braking time)\n")
    w("Tt     = ", Tt,  "  s   (actual braking time)\n")
    w("Ma     = ", Ma,  "  (Alfven Mach number)\n")
    w("V_A0   = ", V_A0,"  cm/s\n")
    w("V_cfl_max    = ", V_cfl_max,    "\n")
    w("V_domain_max = ", V_domain_max, "\n")
    w("Delta  = ", Delta, "  (magnetolaminar parameter)\n")
    w("Eps    = ", Eps,   "  (epsilon = RL/Rb)\n")
    w("\n")

    -- Simulation duration and particle counts
    w("full_time_steps  = ", full_time_steps, "\n")
    w("time_steps       = ", pic_parameters.time_steps, "  (may be capped for PSTD test)\n")
    w("backgr_parts_num = ", backgr_parts_num, "\n")
    w("total_parts_num  = ", total_parts_num, "\n")

    f:flush()
end

return P
