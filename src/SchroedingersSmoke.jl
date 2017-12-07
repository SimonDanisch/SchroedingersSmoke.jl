__precompile__(true)
module SchroedingersSmoke

using GPUArrays
using GPUArrays: gpu_rand

include("MatlabPort.jl")
include("parallel.jl")
include("particles.jl")

end
