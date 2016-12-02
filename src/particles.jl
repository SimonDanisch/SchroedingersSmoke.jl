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
