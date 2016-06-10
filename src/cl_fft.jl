import OpenCL
const cl = OpenCL

import CLFFT
const clfft = CLFFT


const plan_dict = Dict{NTuple{3, Int}, clfft.Plan}()

function getplan(A)
    get!(plan_dict, size(A)) do
        p = clfft.Plan(Complex64, ctx, size(A))
        clfft.set_layout!(p, :interleaved, :interleaved)
        clfft.set_result!(p, :inplace)
        clfft.bake!(p, queue)
        p
    end
end

function Base.fft!(A::cl.CLArray)
    clfft.enqueue_transform(getplan(A), :forward, [queue], A.buffer, nothing)
end
function Base.ifft!(A::cl.CLArray)
    clfft.enqueue_transform(getplan(A), :backward, [queue], A.buffer, nothing)
end

# function fftshift(buffer::CLArray)
#     shift(buffer, size(buffer,1)/2, size(buffer, 2)/2)
# end




const device = first(cl.devices(:gpu))
const ctx    = cl.Context(device)
const queue  = cl.CmdQueue(ctx)

const N = 64
X = rand(Complex64, (N, N, N))
bufX = cl.CLArray(queue, X)

fft!(X)
fft!(bufX)

R = cl.to_host(bufX)
X
@show allclose(R, X; rtol=1e-2, atol=1e-3)
maximum(abs(R - X))
Base.gc()
