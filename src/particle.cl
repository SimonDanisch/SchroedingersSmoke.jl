float3 staggered_velocity(
        __global float* velocity,
        const float3 point,
        const float3 d,
        const int3 tsize,
        const int3 resp,
    ){
    float3 p0, pxp, pyp, pzp, velocity_neighbours, pp, pp2, w, w1, w2, winv1, winv2;
    int3 xyz, p, i, op;

    xyz = get_global_id3();
    p   = mod(point, tsize)
    i   = (int3)(floor(p / d))+ 1 // i is now in range = 1:size(velocity)
    ip  = (int3)(mod(i, resp))  + 1 // i+1, while staying in range

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

    w   = p - ((float3)(i)-1.0) .* d;
    w1  = (float3)(w.z, w.z, w.y);
    w2  = (float3)(w.y, w.x, w.x);
    winv1 = 1-w2;
    winv2 = 1-w1;
    return winv2 * (winv1 * p0 + w2 * pp) + w1 * (winv1 * pp2 + w2 * velocity_neighbours;
}

void StaggeredAdvect(
        __global float* velocity,
        __global float* particle,
        const float dt,
        const float3 d,
        const int3 tsize,
        const int3 res,
    ){
    float3 p, k1, k2, k3, k4;

    p = getindex3(particles, xyz, res);

    k1 = StaggeredVelocity(velocity, p, d, tsize, res);

    k2 = p + (k1*dt/2.0);
    k2 = StaggeredVelocity(velocity, k2, d, tsize, res);

    k3 = p + k2 * dt/2.0;
    k3 = StaggeredVelocity(velocity, k3, d, tsize, res);

    k4 = p + k3 * dt;
    k4 = StaggeredVelocity(velocity, k4, d, tsize, res);

    float3 result = p + dt * (k1 + 2.0*k2 + 2.0*k3 + k4);

    setindex3(particles, result, i, res);
}
