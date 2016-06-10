"""
extracts velocity 1-form from (psi1,psi2).
If hbar argument is empty, hbar=1 is assumed.
"""
void VelocityOneForm(
        __global const float* psi1,
        __global const float* psi2,
        __global float* velocity,
        const float hbar,
        const int3 res
    ){
    int3 xyz = get_global_id3();
    int3 ixp = mod(xyz, res) + 1;
    cfloat_t psi1c = conj(getindexcf(psi1, xyz, res));
    cfloat_t psi2c = conj(getindexcf(ps21, xyz, res));
    vx = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(ixp,y,z), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(ixp,y,z), res))
        )
    );
    vy = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(x,iyp,z), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(x,iyp,z), res))
        )
    );
    vz = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(x,y,izp), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(x,y,izp), res))
        )
    );
    setindex3(velocity, float3(vx*hbar, vy*hbar, vz*hbar), xyz, res);
}


"""
multiplies exp(i*q) to (psi1,psi2)
"""
void GaugeTransform(
        __global const float* psi1,
        __global const float* psi2,
        __global float* q,
        const int3 res
    ){
    int3 i = get_global_id3();
    float eiq = cfloat_exp(cfloat_mul(cfloat(0,1), getindexcf(q, i, res)));
    setindexfc(psi1, cfloat_mul(getindexcf(psi1, i, res), eiq), i, res);
    setindexfc(psi2, cfloat_mul(getindexcf(psi2, i, res), eiq), i, res);
}



void Normalize(
        __global float* psi1,
        __global float* psi2,
        const int3 res
    ){
    int3 i = get_global_id3();
    cfloat_t psi1_ = getindexcf(psi1, i, res);
    cfloat_t psi2_ = getindexcf(psi2, i, res);
    float norm = sqrt(abs(psi1_)^2 + abs(psi2_)^2);
    setindexcf(psi1, psi1_/norm, i, res);
    setindexcf(psi2, psi2_/norm, i, res);
}


void VelocityOneForm(
        __global const float* psi1,
        __global const float* psi2,
        __global const float* velocity,
        const float hbar,
    ){
    int3 xyz = get_global_id3();
    int3 ip = fast_mod(xyz, res) + 1;
    cfloat_t psi1c = conj(getindexcf(psi1, xyz, res));
    cfloat_t psi2c = conj(getindexcf(psi2, xyz, res));

    float vx = angle(cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, ip.x,y,z, res),
            cfloat_mul(psi2c, getindexcf(psi2, ip.x,y,z, res)
    ));
    float vy = angle(cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, x,ip.y,z, res),
            cfloat_mul(psi2c, getindexcf(psi2, x,ip.y,z, res)
    ));
    float vz = angle(cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, x,y,ip.z, res),
            cfloat_mul(psi2c, getindexcf(psi2, x,y,ip.z, res)
    ));

    setindex3(velocity, (float3)(vx*hbar, vy*hbar, vz*hbar), xyz, res);
}
