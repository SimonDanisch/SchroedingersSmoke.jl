


dims = (64, 32, 32)
grid_size = Vec(4, 2, 2)
grid_res = Vec(dims)
hbar = 0.1f0; dt = 1f0/70f0;
d = Vec3f0(grid_size ./ grid_res)

jet_velocity = Vec3f0(1, 0, 0)
nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
nozzle_len = 0.5;
nozzle_rad = 0.5;
n_particles = 10



psi = [(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]


Normalize!(psi)

# constrain velocity
kvec = jet_velocity ./ hbar
omega = sum(jet_velocity.^2) / (2*hbar)


# constrain velocity
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;
amp1 = abs(psi1)
amp2 = abs(psi2)
psi1[isJet] = amp1[isJet] .* exp(1.0im*phase[isJet])
psi2[isJet] = amp2[isJet] .* exp(1.0im*phase[isJet])


for iter = 1:10
    map_idx!(psi, ()) do xyz, psi, _
        x, y, z = xyz
        _psi = psi[x, y, z]
        pos = getindex.(ranges, xyz)
        if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
            phase = sum(kvec .* pos)
            return map(_psi) do p
                abs(p) * exp(1im * phase)
            end
        end
        _psi
    end
    pressure_project!(velocity, psi, idx2, fac, f_tmp, grid_res, d)
end


particle = Particles(
    zeros(Point3f0, 10_000),
    0, Int[]
)

newp = Point3f0[begin
    rt = rand()*2*pi
    Point3f0(
        nozzle_cen[1] - 0.1,
        nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
        nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
    )
end for i=1:n_particles]

append!(particle, newp)


filter!(particle.active) do i
    p = particle.xyz[i]
    all(i-> p[i] in 0:grid_size[i], 1:3)
end

particle_vis = visualize(
    (Circle(Point2f0(0), 0.005f0), particle.xyz),
    indices = particle.active
)
_view(particle_vis, camera=:perspective)

velocity_vis = visualize(
    velocity,
    ranges = Tuple(ranges)
)
#_view(velocity_vis, camera=:perspective)
function in_grid(i, particle = particle)
    p = particle.xyz[i]
    for i=1:3
        p[1] >= 0 && p[1] <= grid_size[i] || return false
    end
    true
end


for i=1:1000
    t = i*dt
    # incompressable schroedinger flow
    schroedinger_flow!(psi, mask)
    Normalize!(psi)
    pressure_project!(velocity, psi, idx2, fac, f_tmp, grid_res, d)

    # constrain velocity
    map_idx!(psi, ()) do xyz, psi, _
        x, y, z = xyz
        _psi = psi[x, y, z]
        pos = getindex.(ranges, xyz)
        if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
            phase = sum(kvec .* pos) - omega*t
            return map(_psi) do p
                abs(p) * exp(1im * phase)
            end
        end
        _psi
    end
    pressure_project!(velocity, psi, idx2, fac, f_tmp, grid_res, d)

    set_arg!(velocity_vis, :rotation, vec(velocity))

    newp = Point3f0[begin
        rt = rand()*2*pi
        Point3f0(
            nozzle_cen[1],
            nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
            nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
        )
    end for i=1:n_particles]

    append!(particle, newp)
    # advect and show particles
    velocity_one_form!(velocity, idx2, grid_res, psi)
    velocity .= (.*).(velocity, Scalar(inv.(d)))
    filter!(in_grid, particle.active)
    particle.length = length(particle.active)
    staggered_advect!(view(particle.xyz, particle.active), velocity, dt, d, grid_res, grid_size)
    set_arg!(particle_vis, :position, particle.xyz)
    set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active))
    yield()
end
