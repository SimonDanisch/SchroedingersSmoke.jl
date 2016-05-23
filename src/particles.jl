"""
 Particles - class of particle which can be advected by staggered velocity
 field on a TorusDEC grid, using RK4 method.
 Velocities are trilinearly interpolated.
"""
type Particles
    #keeplist::Vector{UInt32}
    xyz::Vector{Point3f0}
end

"""
 advect particle positions using RK4 in a grid torus with
 staggered velocity vx,vy,vz, for dt period of time
"""
function StaggeredAdvect(particle, torus, velocity, dt)
    k1 = StaggeredVelocity(
        particle.xyz,
        torus, velocity
    )
    k2 = StaggeredVelocity(
        particle.xyz+k1*dt/2,
        torus,velocity
    )
    k3 = StaggeredVelocity(
        particle.xyz+k2*dt/2,
        torus,velocity
    )
    k4 = StaggeredVelocity(
        particle.xyz+k3*dt,
        torus,velocity
    )
    for (i,p) in enumerate(particle.xyz)
        @inbounds particle.xyz[i] = p + dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
    end
end

"""
for removing particles
"""
function Keep(particle, ind)
    particle.xyz = particle.xyz[ind]
end



@inline function staggered_velocity(point, dinv, d, tsize, resp, velocity)
    p  = mod(point, tsize)
    i  = Point{3, Int}(floor(p.*dinv))+ 1 # i is now in range = 1:size(velocity)
    ip = Point{3, Int}(mod(i, resp))  + 1 # i+1, while staying in range

    p0  = velocity[i[ 1], i[ 2], i[ 3]]

    pxp = velocity[ip[1], i[ 2], i[ 3]]
    pyp = velocity[i[ 1], ip[2], i[ 3]]
    pzp = velocity[i[ 1], i[ 2], ip[3]]

    _xyz = Point3f0(
        velocity[i[ 1], ip[2], ip[3]][1],
        velocity[ip[1], i[ 2], ip[3]][2],
        velocity[ip[1], ip[2], i[ 3]][3]
    )
    pp  = Point3f0(pyp[1], pxp[2], pxp[3])
    pp2 = Point3f0(pzp[1], pzp[2], pyp[3])

    w = p - (Point3f0(i)-1f0).*d
    w1  = Point3f0(w[3], w[3], w[2])
    w2  = Point3f0(w[2], w[1], w[1])
    winv1 = 1-w2
    winv2 = 1-w1

    winv2 .* (winv1 .* p0  + w2 .* pp) + w1 .* (winv1 .* pp2 + w2 .* _xyz)
end
"""
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
function StaggeredVelocity(pf, torus, velocity)
    d = Point3f0(torus.dx, torus.dy, torus.dz)
    dinv = 1f0/d
    ts = Point3f0(torus.sizex, torus.sizey, torus.sizez)
    resp = Point3f0(torus.resx,torus.resy,torus.resz)
    u = similar(pf)
    @simd for i=eachindex(pf)
        @inbounds u[i] = staggered_velocity(pf[i],dinv, d, ts, resp, velocity)
    end
    u
end
