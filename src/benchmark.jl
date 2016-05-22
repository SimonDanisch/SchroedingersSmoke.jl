import OpenCL
const cl = OpenCL
device, ctx, queue = cl.create_compute_context()
using GeometryTypes

function threaded_map(f, v)

end
function div_inner(d_square, res, xyz, velocity, f)
    im  = mod(xyz-2, res) + 1;
    x,y,z = xyz
    _x = getindex(velocity, im[1], y, z)[1];
    _y = getindex(velocity, x, im[2], z)[2];
    _z = getindex(velocity, x, y, im[3])[3];
    v = velocity[x,y,z]
    ff =  (v[1] - _x)/ d_square[1];
    ff += (v[2] - _y)/ d_square[2];
    ff += (v[3] - _z)/ d_square[3];
    f[x,y,z] = ff
    nothing
end
function threaded_div(ds, v, f)
    n = Threads.nthreads()
    @assert length(v) % n == 0
    segment = div(length(v), n)
    res = size(v)
    Threads.@threads for chunk in 1:segment:length(v)
        @simd for i=chunk:(chunk+segment-1)
            @inbounds div_inner(ds, Vec{3,Int}(res...), Vec{3,Int}(ind2sub(res, i)), v, f)
        end
    end
    f
end
function devec_div(ds, res, v, f)
    for z=1:size(v,3), y=1:size(v,2), x=1:size(v,1)
        div_inner(ds, Vec{3,Int}(res...), Vec{3,Int}(x,y,z), v, f)
    end
    f
end
cl_kernel = """
int sub2ind4(int4 D, int4 index){
    int4 i0 = index-1;
    return i0.x+i0.y*D.x + i0.z * D.x * D.y + i0.w * D.x *D.y *D.z;
}
int sub2ind3D(int3 D, int3 index){
    int3 i0 = index-1;
    return i0.x+i0.y*D.x + i0.z * D.x * D.y;
}
float3 getindex3(__global const float* p, int3 xyz, int3 size){
    return vload3(sub2ind3D(size, xyz), p);
}
void setindex3(__global float* p, float3 value, int3 xyz, int3 size){
    vstore3(value, sub2ind3D(size, xyz), p);
}
float getindex(__global const float* p, int3 xyz, int3 size){
    return p[sub2ind3D(size, xyz)];
}
void setindex(__global float* p, float value, int3 xyz, int3 size){
    p[sub2ind3D(size, xyz)] = value;
}

int3 mod(int3 a, int3 m){
    return a - m*convert_int3(floor(convert_float3(a)/convert_float3(m)));
}
__kernel void Div(
        const int3 res,
        const float3 d,
        __global const float* velocity,
        __global float* f
    ){

    int x = get_global_id(0) + 1;
    int y = get_global_id(1) + 1;
    int z = get_global_id(2) + 1;
    int3 xyz = (int3)(x,y,z);
    int3 im  = mod(xyz-2, res) + 1;

    float _x = velocity[sub2ind3D(res, (int3)(im.x, y, z))  ];
    float _y = velocity[sub2ind3D(res, (int3)(x, im.y, z))+1];
    float _z = velocity[sub2ind3D(res, (int3)(x, y, im.z))+2];

    float3 v = getindex3(velocity, xyz, res);

    float ff = (v.x - _x) / d.x;
    ff += (v.y - _y) / d.y;
    ff += (v.z - _z) / d.z;

    setindex(f, ff, xyz, res);
}
"""
p = cl.Program(ctx, source=cl_kernel) |> cl.build!
div_kernel = cl.Kernel(p, "Div")

function test_cl(vol_size, vol_res)
    res = vol_res
    d = Point3f0(vol_size)./Point3f0(res);
    ds = d.^2
    v = rand(Point3f0, res...)

    f = Array(Float32, res...)
    vf32 = reinterpret(Float32, v, (length(v)*3,))
    res32 = Int32[res...]
    d_square = Float32[ds...]


    a = @elapsed devec_div(d_square, res, v, f)

    v_buffer = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=vf32)
    f_buffer = cl.Buffer(Float32, ctx, :w, length(f))

    cl.set_args!(div_kernel, res32, d_square, v_buffer, f_buffer)
    b = @elapsed begin
        cl.enqueue_kernel(queue, div_kernel, res, nothing)
        cl.finish(queue)
    end
    r = cl.read(queue, f_buffer)
    cl.finish(queue)
    if !isapprox(vec(f), r)
        println("shiit ", mean(vec(f) - r), " ", vol_res)
    end
    a,b
end
println(test_cl((4,4,4), (512,512,256)))
#5.78693, 0.15183

t1 = Float64[]
t2 = Float64[]
for i=5:50:300, j=5:50:300, k=5:50:300
    a,b = test_cl((4,4,4), (i,j,k))
    push!(t1, a)
    push!(t2, b)
end

using Vega
range = vec([i*j*k for i=5:50:300, j=5:50:300, k=5:50:300])
x = [range; range]
y = [t1; t2]
group = [["cpu" for i=t1]; ["gpu" for i in t2]]
show(scatterplot(x = x, y = y, group = group))

function s2i(D, index)
    i0 = index-1;
    return i0[1]+i0[2]*D[1] + i0[3] * D[1] * D[2] + i0[4] * D[1] *D[2] *D[3];
end
function s2i(D::NTuple{3, Int}, index)
    i0 = index-1;
    return i0[1] + i0[2]*D[1] + i0[3] * D[1] * D[2];
end
