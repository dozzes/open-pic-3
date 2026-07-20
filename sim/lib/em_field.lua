-- em_field.lua -- Magnetic and electric field initialization
--
-- FIELD GEOMETRY LIBRARY
-- ----------------------
-- Five initialization functions correspond to five physical scenarios:
--
--   uniform          -- constant B and E everywhere (solar-wind benchmark, simple expansion)
--   dipole_pole      -- cloud expands ABOVE the magnetic pole of a distant dipole
--   dipole_equator   -- cloud expands in the EQUATORIAL PLANE of a distant dipole
--   planet_dipole    -- full 3D planetary dipole; planet body is made absorptive
--   center_magnetic_dipole_field    -- dipole at the center of the box, uniform B0 superimposed
--
-- All dipole functions use the Yee C-grid staggering: each B component is
-- placed at the center of the appropriate cell edge, so Bx, By, Bz are
-- evaluated at slightly different positions within the same cell.
--
-- CGS-Gaussian units throughout. B in Gauss, distances in cm.

local P = {}

if _REQUIREDNAME == nil then
    em_field = P
else
    _G[_REQUIREDNAME] = P
end

require("print_mpi")

local pow = math.pow or function(x, y) return x ^ y end

-- ---------------------------------------------------------------------------
-- uniform -- constant B and E everywhere
--
-- Use this for:
--   - cloud expansion into a uniform ambient field (most standard cases)
--   - background flow cases where E = -(v x B)/c gives the motional field
--
-- Parameters: Ex,Ey,Ez [esu/cm], Bx,By,Bz [Gauss]
-- Returns: |B| at center (= sqrt(Bx^2+By^2+Bz^2), same everywhere)
-- ---------------------------------------------------------------------------
function P.uniform(grid, Ex, Ey, Ez, Bx, By, Bz)
    for kx = 0, (grid:size_x() - 1) do
    for ky = 0, (grid:size_y() - 1) do
    for kz = 0, (grid:size_z() - 1) do
        grid:at(kx, ky, kz).E.x = Ex
        grid:at(kx, ky, kz).E.y = Ey
        grid:at(kx, ky, kz).E.z = Ez
        grid:at(kx, ky, kz).B.x = Bx
        grid:at(kx, ky, kz).B.y = By
        grid:at(kx, ky, kz).B.z = Bz
    end
    end
    end
    return math.sqrt(Bx*Bx + By*By + Bz*Bz)
end

-- ---------------------------------------------------------------------------
-- dipole_pole -- cloud above the magnetic POLE of a distant dipole
--
-- Physical picture: the dipole source is placed below the simulation box
-- at z = -dipole_z. Field lines emerge from the south pole at the bottom of
-- the box, converge toward the center at the top, and look like a fan of
-- nearly-parallel z-directed lines at the box center (the cloud location).
--
-- The field magnitude along the z-axis goes as 1/z^3 from the dipole.
-- Bratio = B_bottom / B_top along z, which determines how far below the
-- box the dipole is placed.
--
-- Parameters:
--   Bcloud  [G]  -- B magnitude at the cloud center (= box center along z)
--   Bratio  [-]  -- Bmax/Bmin ratio along the z-axis; determines dipole distance
--                   Typical value: 100-200 (field at bottom is 100-200x stronger)
-- Returns: maximum |B| in the grid
-- ---------------------------------------------------------------------------
function P.dipole_pole(grid, Bcloud, Bratio)
    local h = grid.step

    local center_x = h * math.floor(0.5 * grid:size_x())
    local center_y = h * math.floor(0.5 * grid:size_y())
    local center_z = h * math.floor(0.5 * grid:size_z())

    -- dipole_z: distance from bottom of box to the dipole source.
    -- Derived from Bratio: B(0)/B(Lz) = ((dipole_z+Lz)/dipole_z)^3 = Bratio.
    local dipole_z = grid:length_z() / (pow(Bratio, 1.0/3.0) - 1.0)
    -- Dipole moment p such that Bz(center) = Bcloud on the dipole axis.
    local p = 0.5 * Bcloud * pow(dipole_z + center_z, 3.0)

    local B_max = 0.0

    for kx = 0, (grid:size_x() - 1) do
    for ky = 0, (grid:size_y() - 1) do
    for kz = 0, (grid:size_z() - 1) do
        grid:at(kx, ky, kz).E.x = 0.0
        grid:at(kx, ky, kz).E.y = 0.0
        grid:at(kx, ky, kz).E.z = 0.0

        -- Staggered position for each B component.
        -- Dipole at (center_x, center_y, -dipole_z) in grid coordinates.
        local x = kx*h - center_x
        local y = ky*h - center_y
        local z = kz*h + dipole_z     -- shift so dipole is at z=0 in local coords

        local r2, r5, Bx, By, Bz

        r2 = (x + 0.5*h)*(x + 0.5*h) + y*y + z*z
        r5 = pow(r2, 2.5)
        Bx = 3.0 * (x + 0.5*h) * z * p / r5
        grid:at(kx, ky, kz).B.x = Bx

        r2 = x*x + (y + 0.5*h)*(y + 0.5*h) + z*z
        r5 = pow(r2, 2.5)
        By = 3.0 * (y + 0.5*h) * z * p / r5
        grid:at(kx, ky, kz).B.y = By

        r2 = x*x + y*y + (z + 0.5*h)*(z + 0.5*h)
        r5 = pow(r2, 2.5)
        Bz = (3.0*(z + 0.5*h)*(z + 0.5*h) - r2) * p / r5
        grid:at(kx, ky, kz).B.z = Bz

        B_max = math.max(B_max, math.sqrt(Bx*Bx + By*By + Bz*Bz))
    end
    end
    end

    print_mpi.print_root(proc_idx, "dipole_pole: dipole_z = ", dipole_z,
                         "  B_max = ", B_max, "\n")
    return B_max
end

-- ---------------------------------------------------------------------------
-- dipole_equator -- cloud in the EQUATORIAL PLANE of a distant dipole
--
-- Physical picture: the dipole source is placed to the LEFT of the simulation
-- box (at x = dipole_x < 0). In the equatorial plane, the field is
-- purely along z (Bx=By=0), falling off as 1/r^3 from the dipole.
-- The cloud is at the box center; the dipole is to its left.
--
-- Parameters:
--   Bcloud    [G]  -- B magnitude at the cloud center
--   Bratio    [-]  -- Bmax/Bmin along the x-axis across the box
--   dipole_x  [cm] -- x-coordinate of the dipole (0 = auto-computed from Bratio)
--                    When 0, dipole is placed at distance Lx/(Bratio^(1/3)-1)
--                    to the left of the box.
-- Returns: maximum |B| in the grid
-- ---------------------------------------------------------------------------
function P.dipole_equator(grid, Bcloud, Bratio, dipole_x)
    local h = grid.step

    local center_y = h * math.floor(0.5 * grid:size_y())
    local center_z = h * math.floor(0.5 * grid:size_z())

    local Lx = grid:length_x()
    local B_ratio_m3 = pow(Bratio, 1.0/3.0)

    -- Dipole moment p from the equatorial field at the cloud center.
    local p = -0.125 * Bcloud * pow(Lx * (B_ratio_m3 + 1) / (B_ratio_m3 - 1), 3.0)

    -- Dipole x-position: auto (leftward of box) when dipole_x == 0.
    dipole_x = dipole_x or 0
    if dipole_x == 0 then
        dipole_x = -Lx / (B_ratio_m3 - 1.0)
    end

    local B_max = -math.huge

    for kx = 0, (grid:size_x() - 1) do
    for ky = 0, (grid:size_y() - 1) do
    for kz = 0, (grid:size_z() - 1) do
        grid:at(kx, ky, kz).E.x = 0.0
        grid:at(kx, ky, kz).E.y = 0.0
        grid:at(kx, ky, kz).E.z = 0.0

        local x = kx*h - dipole_x
        local y = ky*h - center_y
        local z = kz*h - center_z

        local r2, r5, Bx, By, Bz

        r2 = (x + 0.5*h)*(x + 0.5*h) + y*y + z*z
        r5 = pow(r2, 2.5)
        Bx = 3.0 * (x + 0.5*h) * z * p / r5
        grid:at(kx, ky, kz).B.x = Bx

        r2 = x*x + (y + 0.5*h)*(y + 0.5*h) + z*z
        r5 = pow(r2, 2.5)
        By = 3.0 * (y + 0.5*h) * z * p / r5
        grid:at(kx, ky, kz).B.y = By

        r2 = x*x + y*y + (z + 0.5*h)*(z + 0.5*h)
        r5 = pow(r2, 2.5)
        Bz = (3.0*(z + 0.5*h)*(z + 0.5*h) - r2) * p / r5
        grid:at(kx, ky, kz).B.z = Bz

        B_max = math.max(B_max, math.sqrt(Bx*Bx + By*By + Bz*Bz))
    end
    end
    end

    print_mpi.print_root(proc_idx, "dipole_equator: dipole_x = ", dipole_x,
                         "  B_max = ", B_max, "\n")
    return B_max
end

-- ---------------------------------------------------------------------------
-- planet_dipole -- full 3D planetary dipole with an absorptive planet body
--
-- Physical picture: a planet with a dipole magnetic field. The planet body
-- (radius R, centered at center_x/y/z) is set to cs_absorptive: particles
-- that enter the body are removed from the simulation. The field is the
-- standard magnetic dipole field (axis along z) outside the body.
--
-- Parameters:
--   B0              [G]   -- field magnitude at the equatorial surface (r=R)
--   R               [cm]  -- planet body radius
--   center_x/y/z    [cm]  -- planet center position
-- Returns: maximum |B| in the grid
-- ---------------------------------------------------------------------------
function P.planet_dipole(grid, B0, R, center_x, center_y, center_z)
    local h  = grid.step
    -- Dipole moment: B_eq(R) = B0, B_eq(r) = B0*(R/r)^3.
    -- On the dipole axis Bz = 2*BR/r^3; at equator Bz = -BR/r^3.
    local BR = B0 * pow(R, 3)

    local B_max = -math.huge

    for i = 0, grid:size_x() - 1 do
    for j = 0, grid:size_y() - 1 do
    for k = 0, grid:size_z() - 1 do
        local x = i*h - center_x
        local y = j*h - center_y
        local z = k*h - center_z

        local r2, r5, Bx, By, Bz

        r2 = (x + 0.5*h)*(x + 0.5*h) + y*y + z*z
        r5 = pow(r2, 2.5)
        Bx = 3*BR*(x + 0.5*h)*z / r5
        grid:at(i,j,k).B.x = Bx

        r2 = x*x + (y + 0.5*h)*(y + 0.5*h) + z*z
        r5 = pow(r2, 2.5)
        By = 3*BR*(y + 0.5*h)*z / r5
        grid:at(i,j,k).B.y = By

        r2 = x*x + y*y + (z + 0.5*h)*(z + 0.5*h)
        r5 = pow(r2, 2.5)
        Bz = BR*(3*(z + 0.5*h)*(z + 0.5*h) - r2) / r5
        grid:at(i,j,k).B.z = Bz

        -- Cells inside the planet body become absorptive boundaries.
        local r = math.sqrt(x*x + y*y + z*z)
        if r < R then
            grid:at(i,j,k).state = Cell.cs_absorptive
        end

        B_max = math.max(B_max, math.sqrt(Bx*Bx + By*By + Bz*Bz))
    end
    end
    end

    print_mpi.print_root(proc_idx, "planet_dipole: R = ", R,
                         "  B_max = ", B_max, "\n")
    return B_max
end

-- ---------------------------------------------------------------------------
-- center_magnetic_dipole_field -- dipole at box center superimposed on a uniform field B0z
--
-- Physical picture: a current loop at the box center creates the dipole;
-- a uniform external field B0z is added on top. The two fields cancel at
-- some distance (magnetopause) and reinforce at others.
-- The dipole moment p is a free parameter (dimensional: G*cm^3).
--
-- Parameters:
--   p       [G*cm^3] -- dipole moment (sets field strength)
--   B0z     [G]      -- superimposed uniform Bz (may be 0)
--   center_x/y/z [cm] -- dipole center position
-- Returns: maximum |B| in the grid
-- ---------------------------------------------------------------------------
function P.center_magnetic_dipole_field(grid, p, B0z, center_x, center_y, center_z)
    local h     = grid.step
    local B_max = 0.0

    for kx = 0, grid:size_x() - 1 do
    for ky = 0, grid:size_y() - 1 do
    for kz = 0, grid:size_z() - 1 do
        local x = kx*h - center_x
        local y = ky*h - center_y
        local z = kz*h - center_z

        local r2, r5, Bx, By, Bz

        r2 = (x + 0.5*h)*(x + 0.5*h) + y*y + z*z
        r5 = pow(r2, 2.5)
        Bx = 3*p*(x + 0.5*h)*z / r5
        grid:at(kx, ky, kz).B.x = Bx

        r2 = x*x + (y + 0.5*h)*(y + 0.5*h) + z*z
        r5 = pow(r2, 2.5)
        By = 3*p*(y + 0.5*h)*z / r5
        grid:at(kx, ky, kz).B.y = By

        r2 = x*x + y*y + (z + 0.5*h)*(z + 0.5*h)
        r5 = pow(r2, 2.5)
        Bz = p*(3*(z + 0.5*h)*(z + 0.5*h) - r2) / r5 + B0z
        grid:at(kx, ky, kz).B.z = Bz

        B_max = math.max(B_max, math.sqrt(Bx*Bx + By*By + Bz*Bz))
    end
    end
    end

    print_mpi.print_root(proc_idx, "center_magnetic_dipole_field: p = ", p,
                         "  B0z = ", B0z, "  B_max = ", B_max, "\n")
    return B_max
end

return P
