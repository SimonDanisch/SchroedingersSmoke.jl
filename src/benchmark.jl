import OpenCL
using GeometryTypes

const cl = OpenCL
device_cpu = first(cl.devices(:cpu))
ctx_cpu     = cl.Context(device_cpu)
queue_cpu   = cl.CmdQueue(ctx_cpu)
device_gpu = first(cl.devices(:gpu))
ctx_gpu     = cl.Context(device_gpu)
queue_gpu   = cl.CmdQueue(ctx_gpu)

function div_inner(d_square, res, xyz, velocity, f)
    im  = mod(xyz-2, res) + 1;
    x,y,z = xyz
    _x = getindex(velocity, im[1], y, z)[1];
    _y = getindex(velocity, x, im[2], z)[2];
    _z = getindex(velocity, x, y, im[3])[3];
    v = velocity[x,y,z]
    ff =  (v - Vec{3, Float32}((_x, _y, _z))) .* d_square;
    f[x,y,z] = sum(ff)
    nothing
end

function threaded_div(ds, res, v, f)
    n = Threads.nthreads()
    @assert length(v) % n == 0
    segment = div(length(v), n)
    res = size(v)
    Threads.@threads for chunk in 1:segment:length(v)
        for i=chunk:(chunk+segment-1)
            @inbounds div_inner(ds, res, Vec{3,Int32}(ind2sub(res, i)), v, f)
        end
    end
    f
end
function devec_div(ds, res, v, f)
    for z=1:size(v,3), y=1:size(v,2), x=1:size(v,1)
        @inbounds div_inner(ds, res, Vec{3,Int32}((x,y,z)), v, f)
    end
    f
end
cl_kernel = readstring("benchkernel.cl")

p_cpu = cl.Program(ctx_cpu, source=cl_kernel) |> cl.build!
kernel_cpu = cl.Kernel(p_cpu, "Div")

p_gpu = cl.Program(ctx_gpu, source=cl_kernel) |> cl.build!
kernel_gpu = cl.Kernel(p_gpu, "Div")

function callthecl(ctx, queue, kernel, f, vf32, d_square, res32)
    v_buffer = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=vf32)
    f_buffer = cl.Buffer(Float32, ctx, :w, length(f))

    cl.set_args!(kernel, res32, d_square, v_buffer, f_buffer)
    b = @elapsed begin
        cl.enqueue_kernel(queue, kernel, (res32...), nothing)
        cl.finish(queue)
    end
    r = cl.read(queue, f_buffer)
    cl.finish(queue)
    r, b
end
function test_cl(vol_size, vol_res)
    res = vol_res

    vf32 = rand(Float32, 3, res...)

    v = reinterpret(Vec{3, Float32}, vf32, res)

    f1 = Array(Float32, res...)
    f2 = Array(Float32, res...)

    d = Vec{3, Float32}(vol_size)./Vec{3, Float32}(res);
    ds = 1./(d.^2)
    res32 = Vec{3,Int32}((res...,))

    res32a = Int32[res...]
    dsa = Float32[ds...]

    r_gpu, t_gpu = callthecl(ctx_gpu, queue_gpu, kernel_gpu, f1, vf32, dsa, res32a)
    r_cpu, t_cpu = callthecl(ctx_cpu, queue_cpu, kernel_cpu, f1, vf32, dsa, res32a)

    b = @elapsed devec_div(ds, res32, v, f2)

    a = @elapsed threaded_div(ds, res32, v, f1)

    println("jl ", a)
    println("jl th ", b)
    println("cl cpu ", t_cpu)
    println("cl gpu ", t_gpu)
    if !isapprox(vec(f2),r_cpu )
        println("shiit, ", mean(vec(f2) - r_cpu ))
    end

end
println(test_cl((10,15,4), (256,256,256)))

# julia> println(test_cl((10,15,4), (256,256,256)))
# jl 4.639070838
# jl th 1.430945919
# cl cpu 0.069358667
# cl gpu 0.033570734


#5.78693, 0.15183

# t1 = Float64[]
# t2 = Float64[]
# for i=5:50:300, j=5:50:300, k=5:50:300
#     a,b = test_cl((4,4,4), (i,j,k))
#     push!(t1, a)
#     push!(t2, b)
# end
#
# using Vega
# range = vec([i*j*k for i=5:50:300, j=5:50:300, k=5:50:300])
# x = [range; range]
# y = [t1; t2]
# group = [["cpu" for i=t1]; ["gpu" for i in t2]]
# show(scatterplot(x = x, y = y, group = group))
#
# function s2i(D, index)
#     i0 = index-1;
#     return i0[1]+i0[2]*D[1] + i0[3] * D[1] * D[2] + i0[4] * D[1] *D[2] *D[3];
# end
# function s2i(D::NTuple{3, Int}, index)
#     i0 = index-1;
#     return i0[1] + i0[2]*D[1] + i0[3] * D[1] * D[2];
# end
