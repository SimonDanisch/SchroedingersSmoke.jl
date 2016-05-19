"""
 Particles - class of particle which can be advected by staggered velocity
 field on a TorusDEC grid, using RK4 method.
 Velocities are trilinearly interpolated.
"""
type Particles
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32} # array of positions
end

"""
 advect particle positions using RK4 in a grid torus with
 staggered velocity vx,vy,vz, for dt period of time
"""
function StaggeredAdvect(particle, torus, vx, vy, vz, dt)
    k1x,k1y,k1z = StaggeredVelocity(
        particle.x,particle.y,particle.z,
        torus,vx,vy,vz
    )
    k2x,k2y,k2z = StaggeredVelocity(
        particle.x+k1x*dt/2,particle.y+k1y*dt/2,particle.z+k1z*dt/2,
        torus,vx,vy,vz
    )
    k3x,k3y,k3z = StaggeredVelocity(
        particle.x+k2x*dt/2,particle.y+k2y*dt/2,particle.z+k2z*dt/2,
        torus,vx,vy,vz
    )
    k4x,k4y,k4z = StaggeredVelocity(
        particle.x+k3x*dt,particle.y+k3y*dt,particle.z+k3z*dt,
        torus,vx,vy,vz
    )
    particle.x = particle.x + dt/6*(k1x+2*k2x+2*k3x+k4x)
    particle.y = particle.y + dt/6*(k1y+2*k2y+2*k3y+k4y)
    particle.z = particle.z + dt/6*(k1z+2*k2z+2*k3z+k4z)
end
"""
for removing particles
"""
function Keep(particle, ind)
    particle.x = particle.x[ind]
    particle.y = particle.y[ind]
    particle.z = particle.z[ind]
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
function StaggeredVelocity(pxf,pyf,pzf,torus,vx,vy,vz)
    px = mod(pxf, torus.sizex)
    py = mod(pyf, torus.sizey)
    pz = mod(pzf, torus.sizez)

    ix = floor(Int, px/torus.dx) + 1
    iy = floor(Int, py/torus.dy) + 1
    iz = floor(Int, pz/torus.dz) + 1
    ixp = mod(ix,torus.resx)+1
    iyp = mod(iy,torus.resy)+1
    izp = mod(iz,torus.resz)+1

    res = (torus.resx,torus.resy,torus.resz)

    ind0 = sub2ind(res, ix,iy,iz)
    indxp = sub2ind(res, ixp,iy,iz)
    indyp = sub2ind(res, ix,iyp,iz)
    indzp = sub2ind(res, ix,iy,izp)
    indxpyp = sub2ind(res, ixp,iyp,iz)
    indypzp = sub2ind(res, ix,iyp,izp)
    indxpzp = sub2ind(res, ixp,iy,izp)

    wx = px - (ix-1)*torus.dx
    wy = py - (iy-1)*torus.dy
    wz = pz - (iz-1)*torus.dz

    ux = ((1-wz).*((1-wy).*vx[ind0 ]+wy.*vx[indyp  ]) +
            wz .*((1-wy).*vx[indzp]+wy.*vx[indypzp]))

    uy = ((1-wz).*((1-wx).*vy[ind0 ]+wx.*vy[indxp  ]) +
            wz .*((1-wx).*vy[indzp]+wx.*vy[indxpzp]))

    uz = ((1-wy).*((1-wx).*vz[ind0 ]+wx.*vz[indxp  ]) +
            wy .*((1-wx).*vz[indyp]+wx.*vz[indxpyp]))
    ux,uy,uz
end
