module SchroedingersSmoke
using GeometryTypes
include("cl_helper.jl")

# matlab compat functions
"""
logical and
"""
function land(a, b)
    map(_land, a, b)
end
_land(a, b) = a != 0 && b != 0

"""
matlab ngrid
"""
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
function ndgrid{T}(vs::AbstractVector{T}...)
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array(T, sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
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
