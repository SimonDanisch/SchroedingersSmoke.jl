
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
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
function StaggeredVelocity(pf, torus, velocity)
    d = (float3)(torus.dx, torus.dy, torus.dz)
    dinv = 1f0/d
    ts = (float3)(torus.sizex, torus.sizey, torus.sizez)
    resp = (float3)(torus.resx,torus.resy,torus.resz)
    u = similar(pf)
    @simd for i=eachindex(pf)
        @inbounds u[i] = staggered_velocity(pf[i],dinv, d, ts, resp, velocity)
    end
    u
end
