__kernel void Div(
        __global const float* velocity,
        __global float* f,

        const int3 res,
        const float3 d
    ){
    float3 vn, v, ff;
    int3 i, im;

    i    = get_global_id3();
    im   = mod(i-2, res) + 1;

    v    = getindex3(velocity, i, res);
    vn.x = getindex3(velocity, (int3)(im.x, i.yz), res).x;
    vn.y = getindex3(velocity, (int3)(i.x, im.y, i.z), res).y;
    vn.z = getindex3(velocity, (int3)(i.xy, im.z), res).z;

    ff = (v - vn) * d;

    setindex(f, ff.x+ff.y+ff.z, i, res);
}

__kernel void StaggeredSharp(
        __global const float* velocity,
        const int3 res,
        const float3 d_inverse
    ){
    int3 xyz = get_global_id3();
    setindex3(velocity, getindex3(velocity, xyz, res) * d_inverse, xyz, res);
}

__kernel void cl_mul(
        __global float* a,
        __global const float* b
    ){
    int i = get_global_id(0);
    a[i] = a[i] * b[i];
}

__kernel void cl_mul_cf(
        __global float* a,
        __global const float* b
    ){
    int i = get_global_id(0);
    float2 af2 = vload2(i, a);
    float2 bf2 = vload2(i, b);
    cfloat_t ac = cfloat_new(af2.x, af2.y);
    cfloat_t bc = cfloat_new(bf2.x, bf2.y);
    cfloat_t r  = cfloat_mul(ac, bc);
    float2 rf2  = (float2)(r.real, r.imag);
    vstore2(rf2, i, a);
}

// __kernel void Div(
//         __global const float* velocity,
//         __global float* f,
//         const int3 res,
//         const float3 d,
//     ){
//     int3 i, im;
//     float3 ff, v;
//
//     i  = get_global_id3();
//     im = mod(i-2, res) + 1;
//
//     float4 xv = getindex4(velocity, (int3)(im.x, i.yz), res);
//
//     v.x = xv.x;
//     v.y = getindex3(velocity, (int3)(i.x, im.y, i.z), res).y;
//     v.z = getindex3(velocity, (int3)(i.xy, im.z), res).z;
//
//     ff = (xv.yzw - v) * d;
//
//     setindex(f, ff.x+ff.y+ff.z, i, res);
// }
