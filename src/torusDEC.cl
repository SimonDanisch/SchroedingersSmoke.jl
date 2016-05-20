

float3 getindex(const float3 * p, int x, int y, int z, int3 size){
    return getindex(p, int3(x, y, z), size)
}
float3 getindex(const float3 * p, int3 xyz, int3 size){
    return p[xyz.x + xyz.y * size.x + xyz.z * size.x * size.y];
}
void setindex(float3 * p, float3 value, int3 xyz, int3 size){
    p[xyz.x + xyz.y * size.x + xyz.z * size.x * size.y] = value;
}



void DerivativeOfTwoForm(
        Torus obj,
        __global const float3* w,
        __global float* f,
    ){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    int3 xyz = int3(x,y,z);
    ip = mod(xyz, obj.res) + 1
    float _x = getindex(w, ip.x, y, z, obj.res)[1];
    float _y = getindex(w, x, ip.y, z, obj.res)[2];
    float _z = getindex(w, x, y, ip.z, obj.res)[3];
    float3 wi = getindex(w, xyz, obj.res);
    ff =  _x - wi[1];
    ff += _y - wi[2];
    ff += _z - wi[3];
    setindex(f, ff, xyz, obj.res);
}

void Div(
        Torus obj,
        __global const float3* velocity,
        __global float* f,
    ){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    int3 im = mod(x-2, obj.res) + 1;
    float _x = getindex(velocity, im[1], y, z, obj.res)[1];
    float _y = getindex(velocity, x, im[2], z, obj.res)[2];
    float _z = getindex(velocity, x, y, im[3], obj.res)[3];
    float3 v = getindex(velocity, x, y, z, obj.res);
    float ff = (v[1] - _x) / obj.dx^2;
    ff += (v[2] - _y) / obj.dy^2;
    ff += (v[3] - _z) / obj.dz^2;
    setindex(f, ff, ,x,y,z, obj.res);
}
//
// """
//  For a 1-form v compute the corresponding vector field v^sharp by
//  averaging to vertices
// """
//
// function Sharp(obj::TorusDEC, velocity)
//     d = 1f0/Point3f0(obj.dx, obj.dy, obj.dz)
//     u = similar(velocity)
//     @inbounds for z=obj.iz, y=obj.iy, x=obj.ix
//         ixm = mod(x-2, obj.resx) + 1
//         iym = mod(y-2, obj.resy) + 1
//         izm = mod(z-2, obj.resz) + 1
//         x = velocity[ixm, y, z][1]
//         y = velocity[x, iym, z][2]
//         z = velocity[x, y, izm][3]
//         u[x,y,z] = 0.5*Point3f0(x,y,z)+velocity[x,y,z].*d
//     end
//     u
// end
//
// """
//  For a 1-form v compute the corresponding vector field `v^sharp` as
//  a staggered vector field living on edges
// """
// function StaggeredSharp(obj::TorusDEC, velocity)
//     d = 1f0/Point3f0(obj.dx, obj.dy, obj.dz)
//     @inbounds for i in eachindex(velocity)
//         velocity[i] = velocity[i] .* d
//     end
//     velocity
// end
//
// """
// PoissonSolve by Spectral method
// """
// function PoissonSolve(obj, f)
//     fc = fft(f)
//     @inbounds for i=1:length(obj.fac)
//         fc[i] = fc[i] .* obj.fac[i]
//     end
//     ifft!(fc)
//     fc
// end
// """
//  For a function f compute the 1-form df
// """
// void DerivativeOfFunction(
//         Torus obj,
//         __global const float* f,
//         __global float* v,
//     )
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     int z = get_global_id(2);
//     ixp = mod(x, obj.resx) + 1
//     iyp = mod(y, obj.resy) + 1
//     izp = mod(z, obj.resz) + 1
//     ff = f[x,y,z]
//     v[x,y,z] = float3(
//         f[ixp,y,z] - ff,
//         f[x,iyp,z] - ff,
//         f[x,y,izp] - ff,
//     )
// end
//
// """
//  For a 1-form v compute the 2-form dv
// """
// function DerivativeOfOneForm(obj::TorusDEC, velocity)
//     w = similar(velocity)
//     @inbounds for z=obj.iz, y=obj.iy, x=obj.ix
//         ixp = mod(x, obj.resx) + 1
//         iyp = mod(y, obj.resy) + 1
//         izp = mod(z, obj.resz) + 1
//
//         x1 = velocity[x, y, izp][1]
//         y1 = velocity[ixp, y, z][2]
//         z1 = velocity[x, iyp, z][3]
//
//         x2 = velocity[x, iyp, z][1]
//         y2 = velocity[x, y, izp][2]
//         z2 = velocity[ixp, y, z][3]
//
//         v = velocity[x,y,z]
//         w[x,y,z] = Point3f0(
//             v[1] - x1 + x2 - v[1],
//             v[2] - y1 + y2 - v[2],
//             v[3] - z1 + z2 - v[3]
//         )
//     end
//     w
// end
