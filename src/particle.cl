void staggered_velocity(
        __global float* point,
        __global float* velocity,
        float3 dinv,
        float3 d,
        int3 tsize,
        int3 resp,
    ){
        
    int3 p = mod(point, tsize)
    int3 i  = (int3)(floor(p.*dinv))+ 1 // i is now in range = 1:size(velocity)
    int3 ip = (int3)(mod(i, resp))  + 1 // i+1, while staying in range

    p0  = velocity[i.x, i.y, i.z]

    pxp = velocity[ip.x, i.y, i.z]
    pyp = velocity[i.x, ip.y, i.z]
    pzp = velocity[i.x, i.y, ip.z]

    _xyz = (float3)(
        velocity[i.x, ip.y, ip.z].x,
        velocity[ip.x, i.y, ip.z].y,
        velocity[ip.x, ip.y, i.z].z
    )
    pp  = (float3)(pyp.x, pxp.y, pxp.z)
    pp2 = (float3)(pzp.x, pzp.y, pyp.z)

    w   = p - ((float3)(i)-1.0).*d
    w1  = (float3)(w.z, w.z, w.y)
    w2  = (float3)(w.y, w.x, w.x)
    winv1 = 1-w2
    winv2 = 1-w1

    winv2 .* (winv1 .* p0  + w2 .* pp) + w1 .* (winv1 .* pp2 + w2 .* _xyz)
}
