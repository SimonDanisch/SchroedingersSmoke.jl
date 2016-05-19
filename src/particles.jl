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
        torus,velocity
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
    particle.xyz = particle.xyz + dt/6*(k1+2*k2+2*k3+k4)
end
"""
for removing particles
"""
function Keep(particle, ind)
    particle.xyz = particle.xyz[ind]
end

"""
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
function StaggeredVelocity(pf,torus,velocity)
    d = Point3f0(torus.dx, torus.dy, torus.dz)
    ts = Point3f0(torus.sizex, torus.sizey, torus.sizez)
    resp = Point3f0(torus.resx,torus.resy,torus.resz)
    u = similar(pf)
    @inbounds for (index, point) in enumerate(pf)
        p = mod(point, ts)
        i = Point{3, Int}(floor(p./d)) + 1
        ip = Point{3, Int}(mod(i, resp))+1

        p0 = velocity[i[1], i[2], i[3]]
        pxp = velocity[ip[1], i[2],  i[3]]
        pyp = velocity[i[1],  ip[2], i[3]]
        pzp = velocity[i[1],  i[2],  ip[3]]
        w = p - (Point3f0(i)-1.).*d
        ux = (
            (1-w[3])*((1-w[2])*p0[1]+w[2]*pyp[1]) +
            w[3] *((1-w[2])*pzp[1]+w[2]*velocity[i[1], ip[2], ip[3]][1])
        )
        uy = (
            (1-w[3])*((1-w[1])*p0[2]+w[1]*pxp[2]) +
            w[3] *((1-w[1])*pzp[2]+w[1]*velocity[ip[1],  i[2],  ip[3]][2])
        )
        uz = (
            (1-w[2]).*((1-w[1]).*p0[3]+w[1].*pxp[3]) +
            w[2] .*((1-w[1]).*pyp[3]+w[1].*velocity[ip[1], ip[2],  i[3]][3])
        )
        u[index] = Point3f0(ux, uy, uz)
    end
    u
end
