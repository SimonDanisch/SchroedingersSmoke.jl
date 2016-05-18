module SchroedingersSmoke

using ParallelAccelerator

_land(a, b) = a > 0 && b > 0
function land(a, b)
    broadcast(_land, a, b)
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
