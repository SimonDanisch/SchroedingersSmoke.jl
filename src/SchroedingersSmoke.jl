module SchroedingersSmoke

#include("MatlabPort.jl")
include("BroadcastPort.jl")
include("ParallelPort.jl")
include("CUDAPort.jl")
println(CUDAPort)
end
