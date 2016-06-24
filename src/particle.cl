float3 staggered_velocity(
        __global float* velocity,
        const float3 point,
        const float3 d,
        const int3 tsize,
        const int3 res
    ){
    
    float3 p0, pxp, pyp, pzp, velocity_neighbours, pp, pp2, w, w1, w2, winv1, winv2;
    int3 xyz, p, i, op, ip;

    xyz = get_global_id3();
    p   = modfi(point, tsize);
    i   = conv_f3i3(floor(convert_float3(p) / d)) + 1; // i is now in range = 1:size(velocity)
    ip  = mod(i, res) + 1; // i+1, while staying in range

    p0  = getindex3(velocity, (int3)(i.x, i.y, i.z), res);

    pxp = getindex3(velocity, (int3)(ip.x, i.y, i.z), res);
    pyp = getindex3(velocity, (int3)(i.x, ip.y, i.z), res);
    pzp = getindex3(velocity, (int3)(i.x, i.y, ip.z), res);

    velocity_neighbours = (float3)(
        getindex3(velocity, (int3)(i.x, ip.y, ip.z), res).x,
        getindex3(velocity, (int3)(ip.x, i.y, ip.z), res).y,
        getindex3(velocity, (int3)(ip.x, ip.y, i.z), res).z
    );
    pp  = (float3)(pyp.x, pxp.y, pxp.z);
    pp2 = (float3)(pzp.x, pzp.y, pyp.z);

    w   = convert_float3(p) - (convert_float3(i)-1.0f) * d;
    w1  = (float3)(w.z, w.z, w.y);
    w2  = (float3)(w.y, w.x, w.x);
    winv1 = 1.0f-w2;
    winv2 = 1.0f-w1;
    return winv2 * (winv1 * p0 + w2 * pp) + w1 * (winv1 * pp2 + w2) * velocity_neighbours;
}

__kernel void StaggeredAdvect(
        __global float* velocity,
        __global float* particle,
        const float dt,
        const float3 d,
        const int3 tsize,
        const int3 res
    ){
    float3 p, k1, k2, k3, k4;
    int3 xyz = get_global_id3();
    p = getindex3(particle, xyz, res);

    k1 = staggered_velocity(velocity, p, d, tsize, res);

    k2 = p + (k1*dt/2.0f);
    k2 = staggered_velocity(velocity, k2, d, tsize, res);

    k3 = p + k2 * dt/2.0f;
    k3 = staggered_velocity(velocity, k3, d, tsize, res);

    k4 = p + k3 * dt;
    k4 = staggered_velocity(velocity, k4, d, tsize, res);

    float3 result = p + dt * (k1 + 2.0f*k2 + 2.0f*k3 + k4);

    setindex3(particle, result, xyz, res);
}
