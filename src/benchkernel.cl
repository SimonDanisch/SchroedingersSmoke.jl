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

    float _x = getindex3(velocity, (int3)(im.x, y, z), res).x;
    float _y = getindex3(velocity, (int3)(x, im.y, z), res).y;
    float _z = getindex3(velocity, (int3)(x, y, im.z), res).z;

    float3 v = getindex3(velocity, xyz, res);

    float3 ff = (v - (float3)(_x, _y, _z)) * d;

    setindex(f, ff.x+ff.y+ff.z, xyz, res);
}
