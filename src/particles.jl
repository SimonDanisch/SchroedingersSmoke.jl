"""
 Particles - class of particle which can be advected by staggered velocity
 field on a TorusDEC grid, using RK4 method.
 Velocities are trilinearly interpolated.
"""
type Particles
    #keeplist::Vector{UInt32}
    xyz::Vector{Point3f0}
end

"""
 advect particle positions using RK4 in a grid torus with
 staggered velocity vx,vy,vz, for dt period of time
"""
function StaggeredAdvect(particle, torus, velocity, dt)
    k1 = StaggeredVelocity(
        particle.xyz,
        torus,velocity
    )
    k2 = StaggeredVelocity(
        particle.xyz+k1*dt/2,
        torus,velocity
    )
    k3 = StaggeredVelocity(
        particle.xyz+k2*dt/2,
        torus,velocity
    )
    k4 = StaggeredVelocity(
        particle.xyz+k3*dt,
        torus,velocity
    )
    particle.xyz = particle.xyz + dt/6*(k1+2*k2+2*k3+k4)
end
"""
for removing particles
"""
function Keep(particle, ind)
    particle.xyz = particle.xyz[ind]
end

"""
Base doesn't have a vectorized form of sub2ind, so we add it!
"""
function Base.sub2ind{N, T<:Integer}(res::NTuple{N, Int}, A::AbstractArray{T, N}...)
    Int[sub2ind(res, i...) for i in zip(A...)]
end

"""
 evaluates velocity at (px,py,pz) in the grid torus with staggered
 velocity vector field vx,vy,vz
"""
function StaggeredVelocity(pf,torus,velocity)
    d = Point3f0(torus.dx, torus.dy, torus.dz)
    ts = Point3f0(torus.sizex, torus.sizey, torus.sizez)
    res = (torus.resx,torus.resy,torus.resz)
    resp = Point3f0(res)

    i = floor(Int, p/d) + 1
    ip = mod(i, resp)+1

    for (i, pf) in enumerate(velocity)
        p = mod(pf, ts)
        i = floor(p/d) + 1
        ip = mod(i, resp)+1
    end
    ind0 = sub2ind(res, i)
    indxp = sub2ind(res, ixp,iy,iz)
    indyp = sub2ind(res, ix,iyp,iz)
    indzp = sub2ind(res, ix,iy,izp)
    indxpyp = sub2ind(res, ixp,iyp,iz)
    indypzp = sub2ind(res, ix,iyp,izp)
    indxpzp = sub2ind(res, ixp,iy,izp)

    w = p - (i-1)*d

    ux = ((1-wz).*((1-wy).*vx[ind0 ]+wy.*vx[indyp  ]) +
            wz .*((1-wy).*vx[indzp]+wy.*vx[indypzp]))

    uy = ((1-wz).*((1-wx).*vy[ind0 ]+wx.*vy[indxp  ]) +
            wz .*((1-wx).*vy[indzp]+wx.*vy[indxpzp]))

    uz = ((1-wy).*((1-wx).*vz[ind0 ]+wx.*vz[indxp  ]) +
            wy .*((1-wx).*vz[indyp]+wx.*vz[indxpyp]))
    ux,uy,uz
end
