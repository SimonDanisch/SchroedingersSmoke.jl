
function StaggeredAdvect(particle, torus, velocity, dt)
    k1 = copy(particle.xyz)
    StaggeredVelocity(k1, torus, velocity)

    k2 = particle.xyz + (k1*dt/2)
    StaggeredVelocity(k2, torus, velocity)

    k3 = particle.xyz+k2*dt/2
    StaggeredVelocity(k3, torus, velocity)

    k4 = particle.xyz+k3*dt
    StaggeredVelocity(k4, torus, velocity)

    sum_staggering(particles, k1, k2, k3, k4, res, dt ./ 6f0)
end


@cl_kernel(program_particle, StaggeredVelocity,
    particles, velocity, dinv, d, tsize, resp
)

@cl_kernel(program_particle, sum_staggering,
    particles,k1,k2,k3,k4,res,dt
)
"""
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
function StaggeredVelocity(particles, torus, velocity, out)
    dinv = 1f0./d
    StaggeredVelocity(
        particles, velocity,
        dinv, torus.d, torus.size, torus.res
    )
    nothing
end
