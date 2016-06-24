type Particles
    xyz::Vector{Point3f0}
    length::Int
    active::Vector{Int}
end
@cl_kernel(StaggeredAdvect, program_particle,
    velocity, particles, dt, d, tsize, res
)

function StaggeredAdvect(obj::ISF, p::Particles, dt)
    StaggeredAdvect(obj.velocity, p.xyz, dt, obj.d, obj.res)
end