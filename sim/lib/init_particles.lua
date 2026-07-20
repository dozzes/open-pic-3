-- init_particles.lua -- Particle placement
--
-- HOW MACRO-PARTICLES WORK
-- ------------------------
-- We cannot simulate every real ion (cloud has ~10^18 of them).
-- Instead each "macro-particle" represents ni = N_real / N_macro real ions.
-- The macro-particle has the same charge-to-mass ratio as a real ion, so
-- it moves identically under the Lorentz force.
--
-- CLOUD INITIALIZATION
-- --------------------
-- The cloud is a sphere of radius R0 centered at (cloud_x, cloud_y, cloud_z).
-- Macro-particles sit on a regular 3D lattice with spacing cloud_part_dist.
-- We use 8-OCTANT SYMMETRY: place one particle in the +++ octant, then
-- mirror it into all 7 other octants.  This gives exact 8-fold symmetry
-- at t=0 and reduces loop iterations by 8x.
--
-- Each particle's initial velocity is proportional to its distance from center:
--   v = v_max * r_vec / R0   (Hubble-like radial expansion)
--
-- BACKGROUND INITIALIZATION
-- -------------------------
-- Background particles fill the box uniformly, excluding:
--   - The cloud sphere (not double-counting those ions)
--   - A thin boundary strip (particles there would immediately be absorbed)

local P = {}

-- Small random displacement to break lattice regularity.
-- Amplitude is half the particle spacing in each direction.
local function random_offset(dist)
    return dist * (math.random() - 0.5)
end
_G.random_offset = random_offset   -- make global so callbacks.lua can access it

-- ============================================================================
-- CLOUD PARTICLES
-- ============================================================================
function P.init_cloud()
    local jitter_dist = 0.5 * cloud_part_dist
    local to_R        = cloud_R0 + cloud_part_dist
    local pk          = 0
    local W0_num      = 0.0

    -- Soft edge: taper each particle's statistical weight from 1 (r <= R_in)
    -- down to 0 at r = cloud_R0, over a shell of width cloud_soft_edge_width_h*h.
    -- The hard cutoff x^2+y^2+z^2 <= cloud_R0_Sq is a step function in density
    -- that seeds a deterministic m=4 perturbation from the coarse discretization
    -- of the sphere on the Cartesian lattice (Analysis Notes 2026-07-04). Tapering
    -- smooths that step without changing which lattice points are placed, so the
    -- particle count (and the total_parts_num sanity check) is unaffected -- only
    -- the weight is renormalized below to keep the total ion count = cloud_ions_num.
    local soft_edge   = cloud_soft_edge_enabled and true or false
    local edge_width  = soft_edge and (cloud_soft_edge_width_h or 1.0) * h or 0.0
    local R_in        = cloud_R0 - edge_width

    local function taper(r)
        if not soft_edge or r <= R_in then return 1.0 end
        local w = (cloud_R0 - r) / edge_width
        return (w > 0.0) and w or 0.0
    end

    -- Controlled mode seed: modulate weight by 1 + an*cos(mn*phi - phase),
    -- ported from the legacy_code `ball` routine (an, mn params in
    -- bin/rt36.sta) to let us seed a KNOWN m=4 amplitude and calibrate the
    -- growth rate directly, instead of relying on the uncontrolled geometric
    -- seed from the hard sphere cutoff (Analysis Notes 2026-07-04).
    -- mn must be even: cos(mn*phi) is then invariant under every octant
    -- mirror (x,y)->(+-x,+-y) (cos(mn*(pi-phi)) = (-1)^mn*cos(mn*phi)), so
    -- mode_weight is a single value per outer-loop iteration, exactly like
    -- taper(r) above -- no need to give place()'s 8 mirrors different
    -- weights. For odd mn this invariance fails, so it is rejected here.
    local mode_seed = cloud_mode_seed_enabled and true or false
    local mode_an   = mode_seed and (cloud_mode_seed_amplitude or 0.0) or 0.0
    local mode_mn   = mode_seed and (cloud_mode_seed_number or 4) or 0
    local mode_phase = cloud_mode_seed_phase or 0.0
    if mode_seed and (mode_mn % 2) ~= 0 then
        error("cloud_mode_seed_number must be even (got " .. mode_mn ..
              "): odd mode numbers are not invariant under octant mirroring")
    end

    local function mode_weight(x, y)
        if not mode_seed then return 1.0 end
        local phi = math.atan(y, x)
        return 1.0 + mode_an * math.cos(mode_mn * phi - mode_phase)
    end

    -- First pass: sum tapered/mode-seeded weight over the lattice (same
    -- octant-mirrored geometry as the placement loop below) to renormalize
    -- cloud_ni so the total represented ion count stays exactly cloud_ions_num.
    local ni_norm = 1.0
    if soft_edge or mode_seed then
        local weight_sum = 0.0
        for x = 0.0, to_R, cloud_part_dist do
        for y = 0.0, to_R, cloud_part_dist do
        for z = 0.0, to_R, cloud_part_dist do
            local r_sq = x*x + y*y + z*z
            if r_sq <= cloud_R0_Sq then
                local mult = 1
                if z ~= 0.0 then mult = mult + 1 end
                if y ~= 0.0 then mult = mult + 1 end
                if x ~= 0.0 then mult = mult + 1 end
                if x ~= 0.0 and y ~= 0.0 then mult = mult + 1 end
                if x ~= 0.0 and z ~= 0.0 then mult = mult + 1 end
                if y ~= 0.0 and z ~= 0.0 then mult = mult + 1 end
                if x ~= 0.0 and y ~= 0.0 and z ~= 0.0 then mult = mult + 1 end
                weight_sum = weight_sum + mult * taper(math.sqrt(r_sq)) * mode_weight(x, y)
            end
        end
        end
        end
        ni_norm = cloud_parts_num / weight_sum
    end

    pic_particle_groups:create_group("cloud", 1.0, 1.0, Diagnostics.save_grid_values)
    local cloud_group_id = pic_particle_groups:get_group("cloud").id

    -- Helper: place one macro-particle at the mirrored position and velocity.
    -- Defined once (not per lattice point) and passed the raw group id plus
    -- scalar components, so it goes through the fast pic_particles:set()
    -- overload -- no per-call group-name lookup, no DblVector userdata
    -- construction, no closure re-allocation (perf, Analysis Notes 2026-07-11).
    local function place(dx, dy, dz, dvx, dvy, dvz, ni)
        pic_particles:set(pk, cloud_group_id, cloud_x + dx, cloud_y + dy, cloud_z + dz, dvx, dvy, dvz, ni)
        pk = pk + 1
    end

    -- Iterate over the +++ octant only; mirror to all 8 octants inside the loop.
    for x = 0.0, to_R, cloud_part_dist do
    for y = 0.0, to_R, cloud_part_dist do
    for z = 0.0, to_R, cloud_part_dist do
        if (x*x + y*y + z*z) <= cloud_R0_Sq then

            -- Marker particles on the z=0 plane let us track the cloud cross-section.
            if z == 0.0 then
                pic_marker_particles:add("cloud_nominal_z=0", pk)
                if (x == 0.0) or (y == 0.0) or (x == y) then
                    pic_marker_particles:add("cloud_nominal_z=0_mrk", pk)
                end
            end

            -- Radial velocity: v proportional to displacement from center.
            local vx = cloud_v_max * x / cloud_R0
            local vy = cloud_v_max * y / cloud_R0
            local vz = cloud_v_max * z / cloud_R0

            -- Optional jitter: shift position slightly off the lattice.
            -- Skip jitter for particles on the axis (x=0, y=0, z=0) to
            -- preserve the exact symmetry at the coordinate planes.
            local xj, yj, zj = x, y, z
            if cloud_jitter_enabled then
                xj = (x == 0.0) and 0.0 or (x + random_offset(jitter_dist))
                yj = (y == 0.0) and 0.0 or (y + random_offset(jitter_dist))
                zj = (z == 0.0) and 0.0 or (z + random_offset(jitter_dist))
                -- If jitter pushed the particle outside the sphere, revert.
                if (xj*xj + yj*yj + zj*zj) > cloud_R0_Sq then xj, yj, zj = x, y, z end
            end

            -- Per-particle weight: uniform cloud_ni, tapered near the boundary
            -- when soft_edge is enabled, and/or modulated by a controlled
            -- azimuthal mode when mode_seed is enabled (see taper()/
            -- mode_weight()/ni_norm above).
            local ni = cloud_ni * ni_norm * taper(math.sqrt(x*x + y*y + z*z)) * mode_weight(x, y)

            -- Place all 8 octant mirrors (guarded to avoid duplicates on planes).
            place( xj,  yj,  zj,  vx,  vy,  vz, ni)   -- octant +++
            if z ~= 0.0 then place( xj,  yj, -zj,  vx,  vy, -vz, ni) end  -- ++-
            if y ~= 0.0 then
                place( xj, -yj,  zj,  vx, -vy,  vz, ni)  -- +-+
                if z == 0.0 then
                    pic_marker_particles:add("cloud_nominal_z=0", pk-1)
                    if (x == 0.0) or (x == y) then
                        pic_marker_particles:add("cloud_nominal_z=0_mrk", pk-1)
                    end
                end
            end
            if y ~= 0.0 and z ~= 0.0 then place( xj, -yj, -zj,  vx, -vy, -vz, ni) end  -- +--
            if x ~= 0.0 then
                place(-xj,  yj,  zj, -vx,  vy,  vz, ni)  -- -++
                if z == 0.0 then
                    pic_marker_particles:add("cloud_nominal_z=0", pk-1)
                    if (y == 0.0) or (x == y) then
                        pic_marker_particles:add("cloud_nominal_z=0_mrk", pk-1)
                    end
                end
            end
            if x ~= 0.0 and z ~= 0.0 then place(-xj,  yj, -zj, -vx,  vy, -vz, ni) end  -- -+-
            if x ~= 0.0 and y ~= 0.0 then
                place(-xj, -yj,  zj, -vx, -vy,  vz, ni)  -- --+
                if z == 0.0 then
                    pic_marker_particles:add("cloud_nominal_z=0", pk-1)
                    if (x == y) then
                        pic_marker_particles:add("cloud_nominal_z=0_mrk", pk-1)
                    end
                end
            end
            if x ~= 0.0 and y ~= 0.0 and z ~= 0.0 then
                place(-xj, -yj, -zj, -vx, -vy, -vz, ni)  -- ---
            end

            -- Accumulate kinetic energy for the energy log (weighted by ni;
            -- reduces to the old uniform-cloud_ni formula when soft_edge is off).
            W0_num = W0_num + (vx*vx + vy*vy + vz*vz) * ni
        end
    end
    end
    end

    -- Total kinetic energy (factor 4 for the 8 octants, divided by 2 for KE).
    W0_num = 4 * W0_num * cloud_mi
    print_mpi.print_root(proc_idx, "  cloud particles:      ", pk, "\n")
    print_mpi.print_root(proc_idx, "  W0_num:               ", W0_num, "\n")
    return pk, W0_num
end

-- ============================================================================
-- CLOUD PARTICLES -- SPHERICAL-SHELL PACKING (alternative to the cubic
-- lattice used by P.init_cloud())
-- ============================================================================
-- Ported from the ancestor legacy_code (subroutine ball,
-- ~2004-2008), ported per user request after reviewing that file (Analysis
-- Notes 2026-07-04). P.init_cloud()'s cubic lattice truncated by a hard
-- sphere has facets at the scale of the lattice spacing -- the root cause of
-- the m=4 seed this whole investigation is about. `ball` avoids that facet
-- entirely: the sphere is built from concentric EQUAL-VOLUME shells, each
-- shell split into latitude bands of equal AREA, each band split into equal
-- ANGULAR steps around the azimuth, with a randomized per-shell starting
-- azimuth (breaking any residual correlation between shells). Points are
-- placed in antipodal pairs (r, -r) rather than octant-mirrored.
--
-- Opt-in via cloud_shell_packing_enabled; P.init_cloud() (the lattice path)
-- remains the default for every existing case.
--
-- Differences from the Fortran original (both deliberate):
--   - Energy is accumulated per actual placed particle (0.5*m*v^2), since
--     this function iterates over the whole sphere directly -- unlike
--     P.init_cloud()'s 4x-at-the-end trick, which only holds because that
--     loop visits one representative particle per 8-fold octant mirror.
--   - Placement hard-stops the moment cloud_parts_num particles are placed
--     (never overshoots), so pk_end matches the pre-computed cloud_parts_num
--     exactly, same as P.init() requires from the lattice path. The
--     equal-volume-shell construction is only approximate in how many
--     particles land in each shell/band (integer rounding of k1/k2), so
--     without this hard stop the natural place-until-full-sphere endpoint
--     would not exactly equal cloud_parts_num.
--   - No marker particles are registered here (P.init_cloud()'s
--     "cloud_nominal_z=0[_mrk]" markers). Confirmed safe: every case in this
--     series sets pic_parameters.save_particle_diagnostics = false, and
--     MarkerParticles iterates only over whatever names were registered
--     (src/save_particles.cpp), so an empty marker set is a no-op. If a
--     future case needs marker-based per-particle diagnostics together with
--     shell packing, markers would need to be added here first.
function P.init_cloud_shells()
    local r0  = cloud_R0
    local jmc = cloud_parts_num
    local sv  = cloud_v_max / r0
    local v0  = 4.0 * math.pi * r0*r0*r0 / (3.0 * jmc)
    local s8  = (v0 / 8.0) ^ (1.0 / 3.0)

    local random_phase = cloud_shell_random_phase_enabled
    if random_phase == nil then random_phase = true end

    local pk     = 0
    local W0_num = 0.0
    local r1     = r0

    pic_particle_groups:create_group("cloud", 1.0, 1.0, Diagnostics.save_grid_values)
    local cloud_group_id = pic_particle_groups:get_group("cloud").id

    -- Outer loop: one concentric equal-volume shell per iteration, shrinking
    -- inward from r1 = cloud_R0 until cloud_parts_num particles are placed.
    while pk < jmc and r1 > 0.0 do
        -- Solve the shell thickness x by fixed-point iteration (same
        -- iteration as legacy_code lines 438-443): the shell of thickness x
        -- at outer radius r1 should have volume ~ v0 per particle.
        local x = s8
        while true do
            local s  = x / r1
            local x1 = s8 * (((1.0 - s)^2) / (1.0 - 2.0*s + (4.0/3.0)*s*s)) ^ (1.0/3.0)
            local s2 = math.abs(x1 - x)
            x = x1
            if s2 <= 1.0e-6 then break end
        end

        local r2 = r1 - x
        local s0 = 4.0 * math.pi * r2 * r2
        local k1 = 2 * math.floor(0.5 + s0 / (8.0 * x * x))
        if k1 < 2 then break end  -- degenerate final shell: stop rather than loop forever

        -- Refine x using the actual integer particle count k1 for this shell
        -- (legacy_code lines 451-455), so shell volumes track cloud_parts_num
        -- as closely as integer rounding allows.
        local v1  = k1 * r0*r0*r0 / jmc
        local s89 = r1*r1*r1 - v1
        if s89 < 0.0 then s89 = 0.0 end
        x  = (r1 - s89 ^ (1.0/3.0)) / 2.0
        r2 = r1 - x

        -- Solve the nominal (equal-area) latitude band width d for this
        -- shell (legacy_code lines 456-463); used only to size k2 per band
        -- below -- the actual band boundary is re-solved exactly per band.
        local t0 = math.pi / 2.0
        local s9 = 4.0 * math.pi / k1
        local d  = 0.0
        while true do
            local d1 = math.sqrt(s9 + d*(d - 2.0*math.sin(d/2.0)))
            local s2 = math.abs(d1 - d)
            d = d1
            if s2 <= 1.0e-6 then break end
        end
        local d_nom = d
        local s7    = k1 / 2.0

        -- Inner loop: latitude bands from the equator (t0 = pi/2) to the
        -- pole (t0 -> 0) in the northern hemisphere only -- each band's
        -- points are placed as antipodal pairs, which cover the southern
        -- hemisphere automatically.
        local i1 = 0
        while i1 < k1 and pk < jmc do
            local s89b = s7 * (math.cos(t0 - d_nom) - math.cos(t0))
            if s89b < 0.6 then s89b = 0.6 end
            local k2 = math.floor(0.5 + s89b)
            if k2 < 1 then k2 = 1 end

            -- Re-solve the exact band width for this integer k2 so
            -- consecutive bands tile the hemisphere without gaps.
            local s = math.cos(t0) + k2 / s7
            if s > 1.0 then s = 1.0 end
            local d_band = t0 - math.acos(s)
            local s4 = r2 * math.sin(t0 - d_band/2.0)
            local s5 = r2 * math.cos(t0 - d_band/2.0)
            s = 2.0 * math.pi / k2

            -- Randomized (or, if disabled, fixed half-step) starting azimuth
            -- per band -- legacy_code's g05cae(x) call, replacing its own
            -- commented-out fixed s3=-s/2 (the original 2004 authors'
            -- comment marks this randomization as deliberate, to avoid a
            -- coherent lattice-locked artifact -- precisely today's m=4 seed).
            local s3 = random_phase and (-s * math.random()) or (-s / 2.0)

            for _ = 1, k2 do
                if pk >= jmc then break end
                s3 = s3 + s
                local s11 = s4 * math.cos(s3)
                local s22 = s4 * math.sin(s3)

                -- Place the antipodal pair (r, -r) sharing this ring position.
                for _ = 1, 2 do
                    if pk >= jmc then break end
                    local vx = sv * s11
                    local vy = sv * s22
                    local vz = sv * s5
                    pic_particles:set(pk, cloud_group_id, cloud_x + s11, cloud_y + s22, cloud_z + s5, vx, vy, vz, cloud_ni)
                    W0_num = W0_num + 0.5 * (vx*vx + vy*vy + vz*vz) * cloud_ni
                    pk = pk + 1
                    s11, s22, s5 = -s11, -s22, -s5
                end
            end

            i1 = i1 + 2*k2
            t0 = t0 - d_band
        end

        r1 = r2 - x
    end

    W0_num = W0_num * cloud_mi
    print_mpi.print_root(proc_idx, "  cloud particles:      ", pk, "\n")
    print_mpi.print_root(proc_idx, "  W0_num:               ", W0_num, "\n")
    return pk, W0_num
end

-- ============================================================================
-- BACKGROUND PARTICLES
-- ============================================================================
function P.init_background(pk_start)
    -- Fill the half-box (again using 8-octant symmetry around the cloud center).
    -- Exclude a thin strip near the boundary and the cloud sphere itself.
    local jitter_dist = 0.5 * backgr_part_dist
    local to_x = 0.5*length_x - backgr_particle_boundary_width_x
    local to_y = 0.5*length_y - backgr_particle_boundary_width_y
    local to_z = 0.5*length_z - backgr_particle_boundary_width_z
    local pk   = pk_start

    pic_particle_groups:create_group("backgr", 1.0, 1.0, Diagnostics.save_grid_values)
    -- FlowBackground only selects the *field* boundary conditions in
    -- bound_cond.lua (forces NP/UP in a boundary shell); it does not touch
    -- particles. Particle replenishment is a separate Lua mechanism
    -- (callbacks.lua: inject_background_particles/add_inflow_plane_x, wired
    -- to on_iteration_begin()): it only injects along X, only when
    -- backgr_vx ~= 0, at a wall-relative coordinate
    -- (backgr_particle_boundary_width_x, NOT relative to cloud_x). It
    -- matches the initial lattice's edge plane only if the cloud sits at
    -- the box's geometric center (cloud_center_shift_x_h == 0); a shifted
    -- cloud combined with backgr_vx ~= 0 would need that coordinate fixed
    -- up. No shipped task currently combines those two, so this is
    -- currently dormant, not exercised.
    pic_particle_groups:set_boundary_kind("backgr", BoundaryKind.FlowBackground)
    local backgr_group_id = pic_particle_groups:get_group("backgr").id
    local backgr_vx, backgr_vy, backgr_vz = backgr_v.x, backgr_v.y, backgr_v.z

    -- Helper: place one background macro-particle (all move with backgr_v).
    -- Defined once outside the lattice loops -- see the analogous cloud
    -- `place` above for why (perf, Analysis Notes 2026-07-11).
    local function place(dx, dy, dz)
        pic_particles:set(pk, backgr_group_id, cloud_x + dx, cloud_y + dy, cloud_z + dz,
                           backgr_vx, backgr_vy, backgr_vz, backgr_ni)
        pk = pk + 1
    end

    for x = 0, to_x, backgr_part_dist do
    for y = 0, to_y, backgr_part_dist do
    for z = 0, to_z, backgr_part_dist do
        -- Skip positions inside the cloud sphere.
        if math.sqrt(x*x + y*y + z*z) > cloud_R0 then

            if z == 0.0 then
                pic_marker_particles:add("backgr_nominal_z=0", pk)
                if x == 0.0 then
                    pic_marker_particles:add("backgr_nominal_z=0_x=0", pk)
                end
            end

            local xj, yj, zj = x, y, z
            if backgr_jitter_enabled then
                xj = (x == 0.0) and 0.0 or (x + random_offset(jitter_dist))
                yj = (y == 0.0) and 0.0 or (y + random_offset(jitter_dist))
                zj = (z == 0.0) and 0.0 or (z + random_offset(jitter_dist))
                -- Revert jitter if the shifted position lands inside the cloud.
                if math.sqrt(xj*xj + yj*yj + zj*zj) <= cloud_R0 then
                    xj, yj, zj = x, y, z
                end
            end

            -- 8 octants (same logic as cloud, but all particles share backgr_v).
            place( xj,  yj,  zj)
            if z ~= 0.0 then place( xj,  yj, -zj) end
            if y ~= 0.0 then place( xj, -yj,  zj) end
            if y ~= 0.0 and z ~= 0.0 then place( xj, -yj, -zj) end
            if x ~= 0.0 then place(-xj,  yj,  zj) end
            if x ~= 0.0 and z ~= 0.0 then place(-xj,  yj, -zj) end
            if x ~= 0.0 and y ~= 0.0 then place(-xj, -yj,  zj) end
            if x ~= 0.0 and y ~= 0.0 and z ~= 0.0 then place(-xj, -yj, -zj) end
        end
    end
    end
    end

    print_mpi.print_root(proc_idx, "  background particles: ", pk, "\n")
    return pk
end

-- ============================================================================
-- PUBLIC ENTRY POINT
-- ============================================================================
function P.init()
    -- Reserve the array for all particles in one allocation.
    pic_particles.size = total_parts_num

    print_mpi.print_root(proc_idx, "\nInitializing particles\n")
    local pk_cloud
    if cloud_shell_packing_enabled then
        pk_cloud = P.init_cloud_shells()
    else
        pk_cloud = P.init_cloud()
    end
    local pk_end = P.init_background(pk_cloud)

    -- Sanity check: particle count must match the pre-computed total.
    if pk_end ~= total_parts_num then
        error("Particle count mismatch: placed " .. pk_end
              .. ", expected " .. total_parts_num)
    end
end

return P
