"""
 Particles - class of particle which can be advected by staggered velocity
 field on a TorusDEC grid, using RK4 method.
 Velocities are trilinearly interpolated.
"""
type Particles
    xyz::Vector{Point3f0}
    length::Int
    active::Vector{Int}
end
function Base.append!(p::Particles, xyz)
    @assert length(xyz) + p.length < length(p.xyz)
    range = (p.length+1):(p.length+length(xyz))
    p.xyz[range] = xyz
    p.length = (p.length+length(xyz))
    append!(p.active, range)
    nothing
end
Base.length(p::Particles) = p.length

"""
 advect particle positions using RK4 in a grid torus with
 staggered velocity vx,vy,vz, for dt period of time
"""
function StaggeredAdvect(particle, torus, velocity, dt)
    d = Point3f0(torus.dx, torus.dy, torus.dz)
    ts = Point3f0(torus.sizex, torus.sizey, torus.sizez)
    res = Vec{3, Int}(torus.resx,torus.resy,torus.resz)
    @inbounds for i in particle.active
        p = particle.xyz[i]
        k1 = staggered_velocity(velocity, p, d, ts, res)

        k2 = p + (k1.*dt/2.0f0);
        k2 = staggered_velocity(velocity, k2, d, ts, res)

        k3 = p + k2 .* dt/2.0f0;
        k3 = staggered_velocity(velocity, k3, d, ts, res)

        k4 = p + k3 .* dt;
        k4 = staggered_velocity(velocity, k4, d, ts, res)

        particle.xyz[i] = p + dt .* (k1 + 2.0f0*k2 + 2.0f0*k3 + k4)
    end
end

"""
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
@inline function staggered_velocity(velocity, point, d, ts, res)
    p   = mod(point, ts)
    i   = Vec{3, Int}(floor(p ./ d) + 1)
    ip  = mod(i, res) + 1

    p0  = velocity[i[1], i[2], i[3]]

    pxp = velocity[ip[1], i[2], i[3]]
    pyp = velocity[i[1], ip[2], i[3]]
    pzp = velocity[i[1], i[2], ip[3]]

    velocity_neighbours = Vec3f0(
        velocity[i[1], ip[2], ip[3]][1],
        velocity[ip[1], i[2], ip[3]][2],
        velocity[ip[1], ip[2], i[3]][3]
    )
    pp  = Vec3f0(pyp[1], pxp[2], pxp[3])
    pp2 = Vec3f0(pzp[1], pzp[2], pyp[3])

    w   = p - (i-1.0f0) .* d
    w1  = Vec3f0(w[3], w[3], w[2])
    w2  = Vec3f0(w[2], w[1], w[1])
    winv1 = 1.0f0-w2;
    winv2 = 1.0f0-w1;
    return Point3f0(winv2 .* (winv1 .* p0 + w2 .* pp) + w1 .* (winv1 .* pp2 + w2) .* velocity_neighbours)
end
