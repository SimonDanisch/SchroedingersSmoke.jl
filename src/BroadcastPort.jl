module BroadcastPort

using GeometryTypes
using StaticArrays

function map_idx!{F, T}(f::F, A::AbstractArray, args::T)
    for x = 1:size(A, 1), y = 1:size(A, 2), z = 1:size(A, 3)
        idx = Vec{3, Int}(x, y, z)
        val = f(idx, A, args)
        @inbounds A[x, y, z] = val
    end
    A
end

@inline function inner_velocity_one_form(i, velocity, res_psi_hbar)
    idx2, res, psi, hbar = res_psi_hbar
    i2 = (idx2[1][i[1]], idx2[2][i[2]], idx2[3][i[3]])
    @inbounds begin
        psi12  = psi[i[1],  i[2],  i[3]]
        psix12 = psi[i2[1], i[2],  i[3]]
        psiy12 = psi[i[1],  i2[2] ,i[3]]
        psiz12 = psi[i[1],  i[2],  i2[3]]
    end
    psi1n = Vec(psix12[1], psiy12[1], psiz12[1])
    psi2n = Vec(psix12[2], psiy12[2], psiz12[2])
    angle.(
        conj(psi12[1]) .* psi1n .+
        conj(psi12[2]) .* psi2n
    ) * hbar
end

function velocity_one_form!(velocity, idx2, res, psi, hbar = 1.0f0)
    map_idx!(inner_velocity_one_form, velocity, (idx2, res, psi, hbar))
end

function gauge_transform!(psi, q)
    broadcast!(psi, psi, q) do psi, q
        eiq = exp(1.0im * (-q))
        (psi[1] * eiq, psi[2] * eiq)
    end
end

function div!(f, res, d, velocity)
    ds = inv.(d .^ 2)
    map_idx!(f, ()) do xyz, f, _
        @inbounds begin
            x, y, z = xyz
            ix, iy, iz = mod.(xyz - 2, res) + 1

            v1 = velocity[x,  y,  z]
            v2 = Vec(
                velocity[ix, y,  z][1],
                velocity[x,  iy, z][2],
                velocity[x,  y, iz][3]
            )
        end
        sum((v1 .- v2) .* ds)
    end
    f
end

function poisson_solve(fac, f)
    fc = fft(f)
    fc .= (*).(fc, fac)
    ifft!(fc)
    fc
end

function pressure_project!(velocity, psi, idx2, fac, f, res, d)
    velocity_one_form!(velocity, idx2, res, psi)
    div!(f, res, d, velocity)
    q = poisson_solve(fac, f)
    gauge_transform!(psi, q)
end

function Normalize!(psi)
    broadcast!(psi, psi) do psi
        a, b = abs(psi[1]), abs(psi[2])
        norm = inv(sqrt(a*a + b*b))
        (psi[1] * norm, psi[2] * norm)
    end
    psi
end
@inline twotuple(a, b) = (a, b)

function schroedinger_flow!(psi, mask)
    # extract single psi values into tmp arrays stored in obj
    psi1 = map(first, psi)
    psi2 = map(last, psi)
    psi1 = fftshift(fft(psi1)); psi2 = fftshift(fft(psi2));
    psi1 .= (*).(psi1, mask)
    psi2 .= (*).(psi2, mask)
    psi1 = ifft!(fftshift(psi1)); psi2 = ifft!(fftshift(psi2));
    psi .= twotuple.(psi1, psi2) # convert back
    psi
end

type Particles
    xyz::Vector{Point3f0}
    length::Int
    active::Vector{Int}
end

function staggered_advect!(particle, velocity, dt, d, res, gridsize)
    map!(particle) do p

        k1 = staggered_velocity(velocity, p, d, gridsize, res)

        k2 = p + k1 .* dt * 0.5f0
        k2 = staggered_velocity(velocity, k2, d, gridsize, res)

        k3 = p + k2 .* dt * 0.5f0
        k3 = staggered_velocity(velocity, k3, d, gridsize, res)

        k4 = p + k3 .* dt
        k4 = staggered_velocity(velocity, k4, d, gridsize, res)

        p .+ dt/6f0 .* (k1 .+ 2f0*k2 .+ 2f0*k3 .+ k4)
    end
end

@inline function staggered_velocity(velocity, point, d, gs, res)
    for i=1:3
        (point[i] > 0 && point[i] < gs[i]) || return point
    end
    p   = (%).(Vec(point), gs)
    i   = Vec{3, Int}(floor.(p ./ d)) + 1
    ip  = (%).(i, res) + 1

    v0  = velocity[i[1], i[2], i[3]]

    pxp = velocity[ip[1], i[2], i[3]]
    pyp = velocity[i[1], ip[2], i[3]]
    pzp = velocity[i[1], i[2], ip[3]]

    vn = Vec3f0(
        velocity[i[1], ip[2], ip[3]][1],
        velocity[ip[1], i[2], ip[3]][2],
        velocity[ip[1], ip[2], i[3]][3]
    )
    pp  = Vec3f0(pyp[1], pxp[2], pxp[3])
    pp2 = Vec3f0(pzp[1], pzp[2], pyp[3])

    w   = p - (i - 1) .* d
    w1  = Vec3f0(w[3], w[3], w[2])
    w2  = Vec3f0(w[2], w[1], w[1])

    return Point3f0(
        (1 - w1) .* ((1 - w2) .* v0 + w2 .* pp) +
        w1 .* ((1 - w2) .* pp2 + w2 .* vn)
    )
end



function Base.append!(p::Particles, xyz)
    @assert length(xyz) + p.length < length(p.xyz)
    range = (p.length + 1):(p.length+length(xyz))
    p.xyz[range] = xyz
    p.length = (p.length+length(xyz))
    append!(p.active, range)
    nothing
end
Base.length(p::Particles) = p.length



end
