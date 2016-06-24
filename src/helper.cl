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

cfloat_t getindexcf(__global const float* p, int3 xyz, int3 size){
    float2 c = vload2(sub2ind3D(size, xyz), p);
    return cfloat_new(c.x, c.y);
}
void setindexcf(__global float* p, cfloat_t value, int3 xyz, int3 size){
    float2 v = (float2)(value.real, value.imag);
    vstore2(v, sub2ind3D(size, xyz), p);
}

int3 mod(int3 a, int3 m){
    return a - m*convert_int3(floor(convert_float3(a)/convert_float3(m)));
}


int3 conv_f3i3(float3 a){
    return (int3)(convert_int(a.x), convert_int(a.y), convert_int(a.z));
}
int3 modfi(float3 a, int3 m){
    return conv_f3i3(a) - m*convert_int3(floor(a/convert_float3(m)));
}

float angle(cfloat_t z){return atan2(z.imag, z.real);}


int3 get_global_id3(){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    return (int3)(x,y,z)+1; // to make porting easier, we use based 1 indexing
}



/*


void cufftShift_3D_slice_kernel(
        __global__ float* data,
        const int N,
        const int zIndex,
        const int threadIdx
    ){
    // 3D Volume & 2D Slice & 1D Line
    int sLine = N;
    int sSlice = N * N;
    int sVolume = N * N * N;

    // Transformations Equations
    int sEq1 = (sVolume + sSlice + sLine) / 2;
    int sEq2 = (sVolume + sSlice - sLine) / 2;
    int sEq3 = (sVolume - sSlice + sLine) / 2;
    int sEq4 = (sVolume - sSlice - sLine) / 2;

    // Thread
    int xThreadIdx = threadIdx.x;
    int yThreadIdx = threadIdx.y;

    // Block Width & Height
    int blockWidth = blockDim.x;
    int blockHeight = blockDim.y;

    // Thread Index 2D
    int xIndex = blockIdx.x * blockWidth + xThreadIdx;
    int yIndex = blockIdx.y * blockHeight + yThreadIdx;

    // Thread Index Converted into 1D Index
    int index = (zIndex * sSlice) + (yIndex * sLine) + xIndex;

    cfloat_t regTemp;

    if (zIndex < N / 2)
    {
        if (xIndex < N / 2)
        {
            if (yIndex < N / 2)
            {
                regTemp = data[index];

                // First Quad
                data[index] = data[index + sEq1];

                // Fourth Quad
                data[index + sEq1] = regTemp;
            }
            else
            {
                regTemp = data[index];

                // Third Quad
                data[index] = data[index + sEq3];

                // Second Quad
                data[index + sEq3] = regTemp;
            }
        }
        else
        {
            if (yIndex < N / 2)
            {
                regTemp = data[index];

                // Second Quad
                data[index] = data[index + sEq2];

                // Third Quad
                data[index + sEq2] = regTemp;
            }
            else
            {
                regTemp = data[index];

                // Fourth Quad
                data[index] = data[index + sEq4];

                // First Quad
                data[index + sEq4] = regTemp;
            }
        }
    }
}

__kernel__ void cufftShift_3D_kernel(float* data, int N, int3 block, int3 grid)
{
    int3 id = get_global_id3();
    cufftShift_3D_slice_kernel(data, N, id.x, id.yz);
}
*/