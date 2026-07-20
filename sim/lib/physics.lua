-- physics.lua -- Plasma physics derivation
--
-- This module takes the raw user parameters from main.lua and derives
-- every quantity needed by the simulation: frequencies, length scales,
-- time step, grid size, and particle counts.
--
-- READING GUIDE FOR STUDENTS
-- --------------------------
-- The key scales in a magnetized plasma are:
--
--   Ion plasma frequency:   Wpi = Z*e * sqrt(4*pi*n / mi)
--   Ion inertial length:    c/Wpi  (how far a wave travels in one oscillation)
--   Ion cyclotron frequency: Wci = Z*e*B / (mi*c)
--   Larmor radius:          RL = v_perp / Wci
--   Alfven speed:           V_A = B / sqrt(4*pi*n*mi)
--
-- The grid step h must resolve all these scales:
--   h <= c/Wpi   (resolve plasma oscillations)
--   h <= RL      (resolve ion orbits)
--
-- The time step tau must satisfy CFL (Courant-Friedrichs-Lewy) stability:
--   tau <= h / V_fast   where V_fast is the fastest wave in the system.

local P = {}

-- Shared scale formulas. Production derivation (P.compute) and the
-- verification tasks (sim/tasks/04_VERIFICATION) must use the same
-- definitions, or verification stops checking what production computes.
function P.alfven_speed(B, dens, mi)
    return B / math.sqrt(4 * math.pi * dens * mi)
end

function P.ion_plasma_frequency(Z, dens, mi)
    local cgs = require("cgs")
    return Z * cgs.e * 2 * math.sqrt(math.pi * dens / mi)
end

function P.ion_cyclotron_frequency(Z, B, mi)
    local cgs = require("cgs")
    return Z * cgs.e * B / (mi * cgs.c)
end

function P.compute(opts)
    -- opts: optional case-specific overrides passed through from
    -- require("sim_core").run(opts). No globals, no _G tables.
    --   opts.h_override: explicitly set h [cm] (controlled resolution scans)
    --   opts.cloud_center_shift_h: {x=..., y=..., z=...}, center shift in cells
    opts = opts or {}
    -- Load constants into module-local variables for clarity.
    local cgs = require("cgs")
    local c   = cgs.c
    local e   = cgs.e
    local me  = cgs.me
    local mp  = cgs.mp

    -- Export constants as globals so other modules (diag, params_log) can use them.
    _G.c   = c
    _G.e   = e
    _G.me  = me
    _G.mp  = mp
    _G.pi  = math.pi

    -- -----------------------------------------------------------------------
    -- ALFVEN SPEED AND BACKGROUND DRIFT
    -- -----------------------------------------------------------------------
    -- The Alfven speed is the fundamental velocity scale for low-frequency
    -- electromagnetic waves in a magnetized plasma (like sound speed but for B).
    --   V_A = B / sqrt(4*pi*n*m)
    B0   = Bz0
    V_A0 = P.alfven_speed(B0, backgr_dens, backgr_mi)

    -- Convert Mach-number drift into physical velocity [cm/s].
    backgr_vx = backgr_flow_Ma_x * V_A0
    backgr_vy = backgr_flow_Ma_y * V_A0
    backgr_vz = backgr_flow_Ma_z * V_A0
    backgr_v  = DblVector(backgr_vx, backgr_vy, backgr_vz)
    backgr_v_abs = math.sqrt(backgr_vx^2 + backgr_vy^2 + backgr_vz^2)

    -- Motional electric field from a flowing background: E = -v x B / c.
    -- For B along Z and flow in XY: Ex = -vy*Bz/c, Ey = vx*Bz/c.
    backgr_Ex0 = -backgr_vy * Bz0 / c
    backgr_Ey0 =  backgr_vx * Bz0 / c
    backgr_Ez0 = 0.0

    -- Cloud expansion velocity (radial, from center outward).
    cloud_v_min = 0.0
    cloud_v_max = Ma * V_A0

    -- -----------------------------------------------------------------------
    -- PLASMA FREQUENCIES AND LENGTH SCALES
    -- -----------------------------------------------------------------------
    -- Background ions:
    --   Wpi = ion plasma frequency (oscillation frequency driven by charge separation)
    backgr_Wpi   = P.ion_plasma_frequency(backgr_Z, backgr_dens, backgr_mi)
    backgr_c_Wpi = c / backgr_Wpi   -- ion inertial length (skin depth)
    backgr_Wci   = P.ion_cyclotron_frequency(backgr_Z, B0, backgr_mi)

    -- Ion Larmor radius: radius of circular orbit in B field.
    RL   = cloud_mi * cloud_v_max * c / (cloud_Z * e * B0)
    RL_b = cloud_v_max / backgr_Wci

    -- Grid step: must resolve both the inertial length and the Larmor radius.
    -- We take the larger of the two constraints (the harder one to resolve).
    -- opts.h_override (absolute cm) or opts.h_over_di (dimensionless) bypasses
    -- this floor for controlled resolution-scan experiments.
    local default_h = math.max(2.0 * backgr_c_Wpi, RL / points_on_RL)
    h = opts.h_override or (opts.h_over_di and opts.h_over_di * backgr_c_Wpi) or default_h
    pic_grid.step = h   -- tell the C++ engine about the chosen grid spacing

    -- Cloud geometry in physical units.
    -- Resolution scans can preserve the default case's physical cloud radius
    -- instead of accidentally shrinking the cloud together with h.
    cloud_R0     = opts.cloud_R0_override
                   or ((opts.preserve_reference_cloud_radius and default_h or h) * h_R0)
    cell_volume  = h * h * h          -- volume of one grid cell [cm^3]
    cloud_volume = (4.0/3.0) * pi * cloud_R0^3
    cloud_dens   = cloud_ions_num / cloud_volume

    -- Cloud plasma frequencies (needed for CFL_tau_loc).
    cloud_Wpi   = 2 * cloud_Z * e * math.sqrt(pi * cloud_dens / cloud_mi)
    cloud_c_Wpi = c / cloud_Wpi
    cloud_Wci   = cloud_Z * e * B0 / (cloud_mi * c)

    -- Maximum cyclotron frequency across all species (sets the shortest timescale).
    Wci_max = math.max(cloud_Wci, backgr_Wci)

    -- -----------------------------------------------------------------------
    -- BOUNDARY WIDTHS
    -- -----------------------------------------------------------------------
    -- Particles near the domain edge are absorbed (cs_absorptive boundary).
    -- We keep 1-2 cells around the edge free of particles so boundary
    -- conditions for the fields have room to interpolate smoothly.
    backgr_flow_boundary_cells_x  = 2
    backgr_flow_boundary_cells_yz = 1
    backgr_flow_boundary_width_x  = backgr_flow_boundary_cells_x  * h
    backgr_flow_boundary_width_y  = backgr_flow_boundary_cells_yz * h
    backgr_flow_boundary_width_z  = backgr_flow_boundary_cells_yz * h

    backgr_particle_boundary_cells_x  = 1
    backgr_particle_boundary_cells_yz = 1
    backgr_particle_boundary_width_x  = backgr_particle_boundary_cells_x  * h
    backgr_particle_boundary_width_y  = backgr_particle_boundary_cells_yz * h
    backgr_particle_boundary_width_z  = backgr_particle_boundary_cells_yz * h

    -- -----------------------------------------------------------------------
    -- ELECTRON THERMODYNAMICS
    -- -----------------------------------------------------------------------
    -- In the hybrid model electrons are a fluid.  Their temperature Te
    -- sets the pressure gradient force on ions and the Spitzer resistivity.
    --
    -- Spitzer resistivity: eta = C_Sp / (kB*Te)^1.5
    -- At very low Te (cold electrons) eta -> 0 and we can skip the full
    -- resistivity calculation entirely (cold_electrons_enabled = true).

    local k_B       = 1.380649e-16   -- Boltzmann constant [erg/K]
    local ln_Lambda = 10.0           -- Coulomb logarithm (typical for laboratory plasma)
    local Te_ref_eV = cold_electrons_enabled and 1.0e-6 or 1.0   -- [eV]

    Spitzer_Te_ref = Te_ref_eV * 11604.5   -- convert eV -> Kelvin

    local C_Spitzer = (4.0/3.0) * math.sqrt(2*pi) * e^2 * ln_Lambda * math.sqrt(me)
    resistivity     = cold_electrons_enabled and 0.0
                      or (C_Spitzer / (k_B * Spitzer_Te_ref)^1.5)

    cloud_Te  = Spitzer_Te_ref   -- initial Te inside cloud
    backgr_Te = Spitzer_Te_ref   -- initial Te in background

    -- Electron thermal velocity (zero for cold electrons).
    V_Pe = cold_electrons_enabled and 0.0
           or math.sqrt(k_B * Spitzer_Te_ref / mp)

    -- -----------------------------------------------------------------------
    -- TIME STEP (CFL CONDITIONS)
    -- -----------------------------------------------------------------------
    -- Several conditions limit tau from above.  We take the minimum and then
    -- apply a safety factor of 0.2.
    --
    --   CFL_tau     : fastest wave must not cross more than one cell per step
    --   PSTD_tau    : spectral Faraday solver stability (Alfven CFL)
    --   Wci_tau     : ion cyclotron period must be resolved (>= 10 steps/cycle)
    --   CFL_tau_loc : local dispersive wave CFL (avoids whistler runaway)
    --   CFL_tau_Pe  : electron thermal speed CFL

    V_max        = math.max(cloud_v_max, V_A0, V_Pe)
    V_cfl_max    = math.max(V_max, backgr_v_abs)
    V_domain_max = math.max(cloud_v_max, backgr_v_abs, V_A0)

    CFL_tau     = 0.5   * h / V_cfl_max          -- wave-speed CFL (FDTD)
    PSTD_tau    = 0.367 * h / V_cfl_max          -- Alfven CFL for PSTD
    Wci_tau     = 0.1   / Wci_max                -- 10 steps per cyclotron period
    CFL_tau_loc = (h / math.max(cloud_c_Wpi, backgr_c_Wpi))^2
                  / (pi * math.sqrt(3) * Wci_max)
    CFL_tau_Pe  = cold_electrons_enabled and math.huge or (0.5 * h / V_Pe)

    -- Choose the binding constraint; PSTD needs the extra spectral CFL.
    if magnetic_field_alg == MagneticFieldAlg.PSTD then
        tau = 0.2 * math.min(CFL_tau, PSTD_tau, Wci_tau, CFL_tau_loc, CFL_tau_Pe)
    else
        tau = 0.2 * math.min(CFL_tau, Wci_tau, CFL_tau_loc, CFL_tau_Pe)
    end

    -- -----------------------------------------------------------------------
    -- CLOUD BRAKING RADII AND TIMES
    -- -----------------------------------------------------------------------
    -- Magnetic braking radius Rb: where cloud kinetic energy equals stored B energy.
    --   (1/2) M V^2 ~ (B^2/8pi) * (4/3 pi Rb^3)  =>  Rb = (6*W0 / B^2)^(1/3)
    -- Gas-dynamic braking radius Rg: where cloud momentum equals displaced background.
    -- The actual braking is at Rt = min(Rb, Rg).

    cloud_mass = cloud_ions_num * cloud_mi
    cloud_W0   = 0.3 * cloud_mass * V_max * V_max   -- effective kinetic energy

    Rb    = (6 * cloud_W0 / (B0 * B0)) ^ (1.0/3.0)
    Rb1   = (3 * cloud_ions_num / (4 * pi * backgr_dens * backgr_Z)) ^ (1.0/3.0)
    Rg    = Rb / Ma ^ (2.0/3.0)
    Rt    = math.min(Rb, Rg)
    Tb    = Rb / V_max   -- magnetic braking time
    Tg    = Rg / V_max   -- gas-dynamic braking time
    Tt    = math.min(Tb, Tg)
    Delta = (Rg / RL) ^ 2   -- magnetolaminar parameter
    Eps   = RL / Rb          -- epsilon (orbit-to-braking radius ratio)

    -- -----------------------------------------------------------------------
    -- PARTICLE SPACING AND COUNTS
    -- -----------------------------------------------------------------------
    -- Macro-particles sit on a regular lattice with spacing h/N_per_step.
    -- Each macro-particle represents ni = n * cell_volume / N_per_step^3 real ions.

    cloud_part_dist  = h / cloud_parts_on_step
    backgr_part_dist = h / backgr_parts_on_step

    require("parts_count")
    cloud_parts_num  = parts_count.get_cloud_parts_num(cloud_R0, cloud_part_dist)
    cloud_ni         = cloud_ions_num / cloud_parts_num
    backgr_ni        = backgr_dens * cell_volume / backgr_parts_on_step^3

    -- -----------------------------------------------------------------------
    -- GRID SIZE
    -- -----------------------------------------------------------------------
    -- The box must be wide enough for the cloud to brake and for the
    -- perturbed region to not yet reach the boundaries. box_break_times
    -- is a separate, optional case global (main.lua) that sizes the box
    -- without changing full_time_steps below -- defaults to break_times,
    -- so every case that doesn't set it is unaffected.
    box_break_times = box_break_times or break_times
    x_nodes = math.floor(2 * V_domain_max * box_break_times * Tt / h + 2)
    if x_nodes % 2 == 0 then x_nodes = x_nodes + 1 end   -- keep odd (symmetric center)
    y_nodes = x_nodes
    z_nodes = x_nodes

    pic_grid:resize(x_nodes, y_nodes, z_nodes)
    pic_grid:set_boundary_state(CellState.cs_absorptive)  -- absorbing walls

    length_x = pic_grid:length_x()
    length_y = pic_grid:length_y()
    length_z = pic_grid:length_z()

    -- Cloud center: geometric center of the grid plus an optional controlled
    -- sub-cell shift. This translates the complete particle configuration
    -- relative to the Cartesian mesh without changing its physical radius,
    -- density, velocity profile, or particle count.
    local center_shift = opts.cloud_center_shift_h or {}
    cloud_center_shift_x_h = tonumber(center_shift.x or center_shift[1] or 0.0)
    cloud_center_shift_y_h = tonumber(center_shift.y or center_shift[2] or 0.0)
    cloud_center_shift_z_h = tonumber(center_shift.z or center_shift[3] or 0.0)
    local center_node = math.floor(0.5 * x_nodes)
    cloud_x = h * (center_node + cloud_center_shift_x_h)
    cloud_y = h * (center_node + cloud_center_shift_y_h)
    cloud_z = h * (center_node + cloud_center_shift_z_h)
    -- Diagnostics still need the nearest integer plane index.
    cloud_x_node = math.floor(cloud_x / h + 0.5)
    cloud_y_node = math.floor(cloud_y / h + 0.5)
    cloud_z_node = math.floor(cloud_z / h + 0.5)
    cloud_R0_Sq  = cloud_R0 * cloud_R0   -- pre-computed for fast inside-sphere tests

    -- Background particle count (full box minus cloud sphere minus boundary strips).
    backgr_parts_num = parts_count.get_backgr_parts_num(
        h, length_x, length_y, length_z,
        cloud_x, cloud_y, cloud_z, cloud_R0,
        backgr_part_dist,
        backgr_particle_boundary_width_x,
        backgr_particle_boundary_width_y,
        backgr_particle_boundary_width_z)

    total_parts_num = cloud_parts_num + backgr_parts_num

    -- -----------------------------------------------------------------------
    -- SIMULATION DURATION
    -- -----------------------------------------------------------------------
    full_time_steps = math.floor(break_times * Tt / tau)

    -- For PSTD quick tests we cap the step count.
    if magnetic_field_alg == MagneticFieldAlg.PSTD and pstd_quick_test_steps ~= nil then
        pic_parameters.time_steps = math.min(full_time_steps, pstd_quick_test_steps)
    else
        pic_parameters.time_steps = full_time_steps
    end

    -- Electron thermal conductivity: set proportional to h^2/tau so it damps
    -- sub-grid temperature fluctuations without suppressing resolved structures.
    if cold_electrons_enabled then
        pic_parameters.electron_thermal_conductivity = 0.0
    else
        pic_parameters.electron_thermal_conductivity = 1e-12 * h * h / tau
    end
end

return P
