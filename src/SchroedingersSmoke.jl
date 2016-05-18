module SchroedingersSmoke

# matlab compat functions
_land(a, b) = a > 0 && b > 0
function land(a, b)
    broadcast(_land, a, b)
end
ndgrid(v::AbstractVector) = copy(v)

function ndgrid(v1::AbstractVector, v2::AbstractVector)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

include("torusDEC.jl")
include("ISF.jl")
include("particles.jl")

export land
export ISF
export TorusDEC
export Normalize
export PressureProject
export Particles
export SchroedingerFlow
export VelocityOneForm
export StaggeredSharp
export StaggeredAdvect
export Keep


end # module
