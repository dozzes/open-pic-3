-- init_field.lua -- Electromagnetic field and electron temperature initialization
--
-- FIELD SHAPES
-- ------------
-- The initial magnetic field geometry is selected by the global variable
-- magnetic_field_shape (set in main.lua, defaults to "uniform"):
--
--   "uniform"        B = (0, 0, Bz0), E from background drift
--   "dipole_pole"    cloud above the magnetic POLE of a distant dipole
--   "dipole_equator" cloud in the EQUATORIAL PLANE of a distant dipole
--   "planet_dipole"  full 3D planetary dipole with absorptive planet body
--   "center_magnetic_dipole_field"  dipole at box center + uniform Bz0 superimposed
--
-- Parameters read from globals (set in main.lua):
--
--   All shapes:
--     Bz0               [G]   -- field magnitude at cloud center
--     backgr_Ex0/Ey0/Ez0      -- motional E field from background drift
--
--   "dipole_pole" and "dipole_equator":
--     dipole_Bratio     [-]   -- Bmax/Bmin across the box (typically 100-200)
--
--   "dipole_equator" only:
--     dipole_x_coord    [cm]  -- dipole x position; nil => auto from Bratio
--
--   "planet_dipole":
--     planet_radius_h   [-]   -- planet radius in units of h (e.g., 10 => R=10*h)
--     planet_at_Rt      [bool]-- if true, planet center at cloud_x + Rt
--                               if false, planet center at cloud center (cloud_x/y/z)
--
--   "center_dipole":
--     dipole_moment     [G*cm^3] -- dipole moment p
--
-- NOTE ON STAGGERED GRID (Yee C-grid)
-- ------------------------------------
-- B components live on cell edges, E on cell faces, Te/eta at cell center.
-- The em_field functions handle the staggering internally (each component
-- is evaluated at its proper sub-cell position).

local P = {}

function P.init()
    local ef = require("em_field")

    local shape = magnetic_field_shape or "uniform"

    -- -----------------------------------------------------------------------
    -- Select and apply the magnetic field initialization.
    -- -----------------------------------------------------------------------
    if shape == "uniform" then
        -- Standard case: uniform B = (0,0,Bz0) + motional E from bg drift.
        ef.uniform(pic_grid, backgr_Ex0, backgr_Ey0, backgr_Ez0,
                             0.0,        0.0,        Bz0)

    elseif shape == "dipole_pole" then
        -- Cloud is above the magnetic pole. Dipole is below the box.
        -- Bratio = B_bottom / B_top along z-axis.
        local Bratio = dipole_Bratio or 100
        ef.dipole_pole(pic_grid, Bz0, Bratio)

    elseif shape == "dipole_equator" then
        -- Cloud is in the equatorial plane. Dipole is to the left of the box.
        -- dipole_x_coord = nil => auto-placed by Bratio.
        local Bratio  = dipole_Bratio   or 100
        local dip_x   = dipole_x_coord  or 0
        ef.dipole_equator(pic_grid, Bz0, Bratio, dip_x)

    elseif shape == "planet_dipole" then
        -- Full planetary dipole with absorptive interior.
        -- Planet body radius = planet_radius_h * h.
        local Rp = (planet_radius_h or 10) * h
        local px, py, pz
        if planet_at_Rt then
            -- Place planet center one braking radius beyond the cloud center.
            px = cloud_x + Rt
        else
            px = cloud_x
        end
        py = cloud_y
        pz = cloud_z
        ef.planet_dipole(pic_grid, Bz0, Rp, px, py, pz)

    elseif shape == "center_magnetic_dipole_field" then
        -- Dipole at box center + uniform Bz0 superimposed.
        local p = dipole_moment or 0.0
        ef.center_magnetic_dipole_field(pic_grid, p, Bz0, cloud_x, cloud_y, cloud_z)

    else
        error("init_field: unknown magnetic_field_shape = '" .. tostring(shape) .. "'")
    end

    -- -----------------------------------------------------------------------
    -- Initialize electron temperature Te and resistivity eta at every cell.
    -- Uses cell-CENTER position (kx+0.5)*h for correct Te staggering.
    -- -----------------------------------------------------------------------
    for kx = 0, (pic_grid:size_x() - 1) do
    for ky = 0, (pic_grid:size_y() - 1) do
    for kz = 0, (pic_grid:size_z() - 1) do
        local x    = (kx + 0.5)*h - cloud_x
        local y    = (ky + 0.5)*h - cloud_y
        local z    = (kz + 0.5)*h - cloud_z
        local cell = pic_grid:at(kx, ky, kz)

        -- Cloud interior: hot electrons; exterior: cool background electrons.
        cell.Te  = (x*x + y*y + z*z < cloud_R0_Sq) and cloud_Te or backgr_Te

        -- Uniform resistivity (zero for cold-electron runs).
        cell.eta = pic_parameters.resistivity
    end
    end
    end
end

return P
