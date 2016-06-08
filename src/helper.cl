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

c_float_t getindexfc(__global const float* p, int3 xyz, int3 size){
    float2 c = vload2(sub2ind3D(size, xyz), p);
    return (c_float_t)(c.x, c.y);
}
void setindexfc(__global float* p, c_float_t value, int3 xyz, int3 size){
    float2 v = (float2)(value.real, value.imag);
    vstore2(v, sub2ind3D(size, xyz), p);
}

int3 mod(int3 a, int3 m){
    return a - m*convert_int3(floor(convert_float3(a)/convert_float3(m)));
}
float angle(cfloat_t z){return atan2(z.imag, z.real);}


int3 get_global_id3(){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    return int3(x,y,z)+1; // to make porting easier, we use based 1 indexing
}
