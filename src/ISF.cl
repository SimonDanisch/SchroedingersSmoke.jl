typedef struct _ISF{
  float hbar;
  float dt;
}ISF;

typedef struct _Torus{
  float hbar;
  float dt;
}Torus;



"""
extracts velocity 1-form from (psi1,psi2).
If hbar argument is empty, hbar=1 is assumed.
"""
void VelocityOneForm(
        Torus obj,
        __global const float* psi1,
        __global const float* psi2,
        __global float* velocity,
        float hbar
    ){
    int3 xyz = get_global_id3();
    int3 ixp = mod(xyz, obj.res) + 1;
    cfloat_t psi1c = conj(getindexcf(psi1, xyz, obj.res));
    cfloat_t psi2c = conj(getindexcf(ps21, xyz, obj.res));
    vx = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(ixp,y,z), obj.res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(ixp,y,z), obj.res))
        )
    );
    vy = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(x,iyp,z), obj.res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(x,iyp,z), obj.res))
        )
    );
    vz = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(x,y,izp), obj.res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(x,y,izp), obj.res))
        )
    );
    setindex3(velocity, float3(vx*hbar, vy*hbar, vz*hbar), xyz, obj.res);
}


"""
multiplies exp(i*q) to (psi1,psi2)
"""
void GaugeTransform(
        const Torus obj,
        __global const float* psi1,
        __global const float* psi2,
        __global float* q,
    ){
    int3 i = get_global_id3();
    float eiq = cfloat_exp(cfloat_mul(cfloat(0,1), getindexcf(q, i, obj.res)));
    setindexfc(psi1, cfloat_mul(getindexcf(psi1, i, obj.res), eiq), i, obj.res);
    setindexfc(psi2, cfloat_mul(getindexcf(psi2, i, obj.res), eiq), i, obj.res);
}



void Normalize(
        __global float* psi1,
        __global float* psi2
    ){
    int3 i = get_global_id3();
    cfloat_t psi1_ = getindexcf(psi1, i, obj.res);
    cfloat_t psi2_ = getindexcf(psi2, i, obj.res);
    float norm = sqrt(abs(psi1_)^2 + abs(psi2_)^2);
    psi1[i] = psi1_/norm;
    psi2[i] = psi2_/norm;
}
