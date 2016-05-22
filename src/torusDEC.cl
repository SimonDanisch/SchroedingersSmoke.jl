

int3 get_global_id3(){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    return int3(x,y,z)+1; // to make porting easier, we use based 1 indexing
}
void DerivativeOfTwoForm(
        Torus obj,
        __global const float3* w,
        __global float* f,
    ){

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

__kernel void Div(
        const int3 res,
        const float3 d,
        __global const float* velocity,
        __global float* f
    ){

    int3 xyz = get_global_id3();
    int3 im  = mod(xyz-2, res) + 1;
    float _x = getindex3(velocity, (int3)(im.x, y, z), res).x;
    float _y = getindex3(velocity, (int3)(x, im.y, z), res).y;
    float _z = getindex3(velocity, (int3)(x, y, im.z), res).z;
    float3 v = getindex3(velocity, xyz, res);

    float ff = (v.x - _x) / d.x;
    ff += (v.y - _y) / d.y;
    ff += (v.z - _z) / d.z;

    setindex(f, ff, xyz, res);
}

__kernel void StaggeredSharp(
        TorusDEC obj,
        const float3 d_inverse,
        __global const float* velocity
    ){
    int3 xyz = get_global_id3();
    setindex3(velocity, getindex3(velocity, xyz, obj.res) .* d, xyz, obj.res);
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
