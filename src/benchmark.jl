using GeometryTypes, BenchmarkTools
function Div1(ds, res, velocity, f)
    @inbounds for z=size(f,3), y=size(f,2), x=size(f,1)
        ixm = mod(x-2, res[1]) + 1
        iym = mod(y-2, res[2]) + 1
        izm = mod(z-2, res[3]) + 1
        _x = velocity[ixm, y, z][1]
        _y = velocity[x, iym, z][2]
        _z = velocity[x, y, izm][3]
        v  = velocity[x,y,z]
        ff =  (v[1] - _x) / ds[1]^2
        ff += (v[2] - _y) / ds[2]^2
        ff += (v[3] - _z) / ds[3]^2
        f[x,y,z] = ff
    end
    f
end
function Div2(ds, res, vx, vy, vz, f)
    @inbounds for z=size(f,3), y=size(f,2), x=size(f,1)
        ixm = mod(x-2, res[1]) + 1
        iym = mod(y-2, res[2]) + 1
        izm = mod(z-2, res[3]) + 1
        _x = vx[ixm, y, z]
        _y = vy[x, iym, z]
        _z = vz[x, y, izm]
        ff =  (vx[x,y,z] - _x) / ds[1]^2
        ff += (vy[x,y,z] - _y) / ds[2]^2
        ff += (vz[x,y,z] - _z) / ds[3]^2
        f[x,y,z] = ff
    end
    f
end
function test(vol_size, vol_res)
    res = vol_res
    d = Point3f0(vol_size)./Point3f0(res);
    ds = d.^2
    f = Array(Float64, (res...))
    v = rand(Point3f0, res...)
    vx = rand(Float32, res...)
    vy = rand(Float32, res...)
    vz = rand(Float32, res...)
    println("1 ", @benchmark Div1($ds, $res, $v, $f))
    println("2 ", @benchmark Div2($ds, $res, $vx, $vy, $vz, $f))
end
cl_kernel = """
float3 getindex(const float3 * p, int x, int y, int z, int3 size){
    return getindex(p, int3(x, y, z), size)
}
float3 getindex(const float3 * p, int3 xyz, int3 size){
    return p[xyz.x + xyz.y * size.x + xyz.z * size.x * size.y];
}
void setindex(float3 * p, float3 value, int3 xyz, int3 size){
    p[xyz.x + xyz.y * size.x + xyz.z * size.x * size.y] = value;
}
__kernel void Div(
        int3 res,
        __global const float3* velocity,
        __global float* f,
        float3 d_square
    ){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    int3 xyz = int3(xyz);
    int3 im = mod(xyz-2, res) + 1;
    float _x = getindex(velocity, im[1], y, z, res)[1];
    float _y = getindex(velocity, x, im[2], z, res)[2];
    float _z = getindex(velocity, x, y, im[3], res)[3];
    float3 v = getindex(velocity, xyz, res);
    float ff = (v[1] - _x) / d_square;
    ff += (v[2] - _y) / d_square;
    ff += (v[3] - _z) / d_square;
    setindex(f, ff, xyz, res);
}
"""
import OpenCL
const cl = OpenCL
device, ctx, queue = cl.create_compute_context()
p = cl.Program(ctx, source=sum_kernel_src) |> cl.build!
sum_kernel = cl.Kernel(p, "Div")

test((4,2,2), (256,256,256))
versioninfo()
