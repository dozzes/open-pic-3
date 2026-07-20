-- diag.lua -- Per-step diagnostics
--
-- WHAT WE MEASURE
-- ---------------
-- At every save step we compute integrated energies and write them to
-- energy.txt.  This lets you monitor whether the simulation is physical
-- (energy should be conserved up to numerical errors) and track how the
-- cloud slows down and transfers energy to the background field.
--
-- We also archive the raw .dat grid files into a zip so the working
-- directory stays clean between saves.
--
-- ENERGY DEFINITIONS (normalized by cloud_W0 = initial cloud kinetic energy)
--   W_mf  = integral of B^2 / (8*pi) over the box  (magnetic field energy)
--   W_e   = integral of n_e * m_e * Ue^2 / 2        (electron kinetic energy)
--   W_cloud and W_backgr columns here are legacy zeros. The particle kinetic
--   energies are written by the C++ side to energy_kin.txt (raw erg, one row
--   per save step) -- accumulated inside the particle push loop at v(n+1/2),
--   serial runs only. Normalize by cloud_W0 (see lua_params.txt) at analysis.

local P = {}
local utils = require("utils")

function P.on_step(t, st)
    -- ------------------------------------------------------------------
    -- Compute field energies by summing over all grid cells.
    -- ------------------------------------------------------------------
    local wmf = 0.0   -- magnetic field energy (unnormalized)
    local wel = 0.0   -- electron kinetic energy (unnormalized)
    local wth = 0.0   -- electron thermal energy (unnormalized)

    for kx = 0, (pic_grid:size_x() - 1) do
    for ky = 0, (pic_grid:size_y() - 1) do
    for kz = 0, (pic_grid:size_z() - 1) do
        local cell = pic_grid:at(kx, ky, kz)
        local b    = cell.B:abs()    -- |B| at this cell
        local ue   = cell.UE:abs()   -- |UE| (electron bulk speed)
        wmf = wmf + b  * b
        wel = wel + ue * ue * cell.NP   -- weight by ion density as proxy for n_e
        wth = wth + cell.Te * cell.NP   -- for (3/2) n k_B Te (quasi-neutral: n_e = NP)
    end
    end
    end

    -- Normalize: W_mf = integral B^2/(8*pi) dV, scaled by cloud_W0.
    -- The integral needs the cell volume h^3: without it the sum is not an
    -- energy. NOTE: energy.txt files written before 2026-07-04 lack this
    -- factor (all values uniformly smaller by h^3) -- ratios between runs
    -- on the same grid were unaffected.
    local cell_vol = h * h * h
    wmf = 0.125 * wmf * cell_vol / pi / cloud_W0
    -- W_e = integral (1/2)*n*me*Ue^2 dV, scaled by cloud_W0.
    wel = 0.5   * me  * wel * cell_vol / cloud_W0
    -- W_the = integral (3/2)*n*k_B*Te dV, scaled by cloud_W0. At beta_e ~ 4
    -- this reservoir is ~30x cloud_W0, and grad(Pe) exchanges energy with the
    -- ions (W_cloud rises early in finite-Te runs) -- the balance
    -- W_cloud + W_backgr + W_mf + W_e closes only together with W_the.
    local k_B = 1.380649e-16   -- [erg/K]; Te on the grid is in Kelvin
    wth = 1.5 * k_B * wth * cell_vol / cloud_W0

    print_mpi.print_root(proc_idx, "  W_mf = ", wmf, "  W_e = ", wel,
                         "  W_the = ", wth, "\n")

    -- Write one row to energy.txt (append mode).
    local f   = io.open("energy.txt", "a")
    local tnd = t * pic_parameters.tau / pic_parameters.T_scale  -- normalized time
    if t == 1 then f:write("t\tW_cloud\tW_backgr\tW_mf\tW_e\tW_the\n") end
    -- W_cloud and W_backgr are placeholders (0) until particle-loop energy
    -- tracking is re-enabled in this function.
    f:write(tnd, "\t", 0.0, "\t", 0.0, "\t", wmf, "\t", wel, "\t", wth, "\n")
    f:flush()

    -- ------------------------------------------------------------------
    -- Archive .dat output files every save interval.
    -- ------------------------------------------------------------------
    if (t % st == 0) or (t == st) then
        -- utils.clear_cmd()  -- commented out: let user control terminal clearing
        io.write("Step " .. t .. ": archiving diagnostic data...\n")
        utils.ensure_dir("./diag")

        -- Copy the planes through the cloud center into the diag/ folder.
        utils.copy_pattern("*markers*.dat",                          "./diag")
        utils.copy_pattern("at*.dat",                               "./diag")
        utils.copy_pattern("*_x_" .. cloud_x_node .. "_grd_*.dat", "./diag")
        utils.copy_pattern("*_y_" .. cloud_y_node .. "_grd_*.dat", "./diag")
        utils.copy_pattern("*_z_" .. cloud_z_node .. "_grd_*.dat", "./diag")

        -- Live Q4 tracking: run the square-mode metric on the central
        -- z-slice just copied into diag/ and append one row to
        -- square_mode_live.csv in the task dir. Lets you watch the m=4
        -- grid mode grow during the run (tail the csv) instead of
        -- post-processing after it finishes. Failure is non-fatal (no
        -- python on PATH, front not yet closed at early steps -- the
        -- script simply appends nothing).
        local tools_dir = P.find_tools_dir()
        if tools_dir then
            -- Step is zero-padded to the width of the total step count,
            -- matching create_out_file_name() in io_utilities.cpp:26
            -- (fmt::format("{:0{}}", step, width) with width computed from
            -- pic_parameters.time_steps) -- NOT a fixed-width %d.
            local width = #tostring(pic_parameters.time_steps)
            local slice = string.format("./diag/all_z_%s_grd_%0" .. width .. "d.dat",
                                        tostring(cloud_z_node), t)
            os.execute(string.format(
                'python "%s/measure_square_mode.py" "%s" --center-x %.17g --center-y %.17g --append square_mode_live.csv',
                tools_dir, slice, cloud_x / pic_parameters.L_scale, cloud_y / pic_parameters.L_scale))
        end

        -- Pack all .dat files into a zip archive and remove the originals
        -- to avoid filling the disk with thousands of files.
        local archive = string.format("7z a -tzip -mx1 -mmt16 %d.zip *.dat", t)
        local ok = os.execute(archive)
        if ok then os.execute(utils.RM .. " *.dat") end
    end
end

-- Locate the repo's tools/ directory from package.path: tasks add
-- "<...>/lib/?.lua", and tools/ sits at lib/../../tools (sim/lib ->
-- repo root -> tools). Returns nil if no lib entry is present, which
-- disables live Q4 tracking without breaking the run.
function P.find_tools_dir()
    -- Tasks nested under 04_VERIFICATION/ carry TWO "*lib/?.lua" entries
    -- in package.path (their own lib/ plus sim/lib/) -- the first match
    -- is not necessarily sim/lib, so lib_dir/../../tools can point
    -- anywhere. Check each candidate for the actual script instead of
    -- trusting entry order.
    for entry in package.path:gmatch("[^;]+") do
        local lib_dir = entry:match("^(.*lib)[/\\]%?%.lua$")
        if lib_dir then
            local candidate = lib_dir .. "/../../tools"
            local f = io.open(candidate .. "/measure_square_mode.py", "r")
            if f then
                f:close()
                return candidate
            end
        end
    end
    return nil
end

return P
