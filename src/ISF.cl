__kernel void VelocityOneForm(
        __global const float* psi1,
        __global const float* psi2,
        __global float* velocity,
        const float hbar,
        const int3 res
    ){
    float vx, vy, vz;
    int3 xyz = get_global_id3();
    int3 ip = mod(xyz, res) + 1;
    cfloat_t psi1c = cfloat_conj(getindexcf(psi1, xyz, res));
    cfloat_t psi2c = cfloat_conj(getindexcf(psi2, xyz, res));
    vx = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(ip.x, xyz.yz), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(ip.x, xyz.yz), res))
        )
    );
    vy = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(xyz.x,ip.y,xyz.z), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(xyz.x,ip.y,xyz.z), res))
        )
    );
    vz = angle(
        cfloat_add(
            cfloat_mul(psi1c, getindexcf(psi1, (int3)(xyz.xy, ip.z), res)),
            cfloat_mul(psi2c, getindexcf(psi2, (int3)(xyz.xy, ip.z), res))
        )
    );
    setindex3(velocity, (float3)(vx*hbar, vy*hbar, vz*hbar), xyz, res);
}

__kernel void GaugeTransform(
        __global float* psi1,
        __global float* psi2,
        __global float* q,
        const int3 res
    ){
    int3 i = get_global_id3();
    cfloat_t eiq = cfloat_exp(cfloat_mul(cfloat_new(0.,1.), getindexcf(q, i, res)));
    setindexcf(psi1, cfloat_mul(getindexcf(psi1, i, res), eiq), i, res);
    setindexcf(psi2, cfloat_mul(getindexcf(psi2, i, res), eiq), i, res);
}


__kernel void Normalize(
        __global float* psi1,
        __global float* psi2,
        const int3 res
    ){
    int3 i = get_global_id3();
    cfloat_t psi1_ = getindexcf(psi1, i, res);
    cfloat_t psi2_ = getindexcf(psi2, i, res);

    float a = cfloat_abs(psi1_); 
    float b = cfloat_abs(psi2_);

    float norm = 1.0/sqrt(a*a + b*b);

    setindexcf(psi1, cfloat_mulr(psi1_, norm), i, res);
    setindexcf(psi2, cfloat_mulr(psi2_, norm), i, res);
}

/*
void VelocityOneForm(
        __global const float* psi1,
        __global const float* psi2,
        __global float* velocity,
        const float hbar
    ){
    int3 xyz = get_global_id3();
    int3 ip = fast_mod(xyz, res) + 1;
    cfloat_t psi1c = cfloat_conj(getindexcf(psi1, xyz, res));
    cfloat_t psi2c = cfloat_conj(getindexcf(psi2, xyz, res));

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
*/