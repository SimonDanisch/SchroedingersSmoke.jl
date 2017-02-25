module BroadcastPort

using StaticArrays
const Vec = SVector
const Point = SVector
const Vec3f0 = Vec{3, Float32}
const Point3f0 = Point{3, Float32}

# lots of parameters. To lazy to write out types,
# still don't want to waste performance
type ISF{IntType, FloatType}
    grid_res::Vec{3, IntType}
    physical_size::Vec{3, IntType}
    d::Vec{3, FloatType}
    velocity::Array{Vec{3, FloatType}, 3}
    # gridpoints in pysical size
    ranges::Vec{3, Vector{FloatType}}


    hbar::FloatType             # reduced Planck constant
    dt::FloatType               # time step
    mask::Array{Complex{FloatType}, 3} # Fourier coefficient for solving Schroedinger eq
    # precalculated values
    idx_shifted::Vec{3, Vector{IntType}} # indices shifted by one
    fac::Array{Complex{FloatType}, 3}
    # temporary allocations for intermediates
    f::Array{FloatType, 3}

    function ISF(physical_size, dims, hbar, dt)

        VT = Vec{3, FloatType}
        grid_res = Vec(dims)
        physical_size = Vec(physical_size)
        d = VT(physical_size ./ grid_res)
        fac = zeros(FloatType, dims)
        mask = zeros(Complex{FloatType}, dims)
        f = -4 * pi^2 * hbar
        for x = 1:dims[1], y = 1:dims[2], z = 1:dims[3]
            xyz = Vec(x, y, z)
            # fac
            s = sin(pi * (xyz - 1) ./ grid_res) ./ d
            denom = sum(s .^ 2)
            fac[x,y,z] = -0.25f0 ./ denom

            # schroedingers mask
            k = (xyz - 1 - grid_res ./ 2) ./ physical_size
            lambda = f * sum(k .^ 2)
            mask[x,y,z] = exp(1.0im * lambda * dt / 2)
        end
        fac[1,1,1] = 0

        f = zeros(FloatType, dims)
        velocity = zeros(VT, dims)

        ranges = Vec(ntuple(3) do i
             ((1:dims[i]) - 1) * d[i]
        end)
        idx_shifted = Vec(ntuple(3) do i
            mod(1:grid_res[i], grid_res[i]) + 1
        end)
        new(
            grid_res,
            physical_size,
            d,
            velocity,
            ranges,
            hbar,
            dt,
            mask,
            idx_shifted,
            fac,
            f
        )
    end
end

# helper function to have a map but with random index manipulations
function map_idx!{F, T}(f::F, A::AbstractArray, args::T)
    for x = 1:size(A, 1), y = 1:size(A, 2), z = 1:size(A, 3)
        idx = Vec{3, Int}(x, y, z)
        val = f(idx, A, args)
        @inbounds A[x, y, z] = val
    end
    A
end

@inline function inner_velocity_one_form(i, velocity, idx_psi_hbar)
    idx2, psi, hbar = idx_psi_hbar
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

function velocity_one_form!(isf, psi, hbar = 1.0f0)
    arg = (isf.idx_shifted, psi, hbar)
    map_idx!(inner_velocity_one_form, isf.velocity, arg)
end


function gauge_transform!(psi, q)
    broadcast!(psi, psi, q) do psi, q
        eiq = exp(1.0im * (-q))
        (psi[1] * eiq, psi[2] * eiq)
    end
end

function div!(isf)
    ds = inv.(isf.d .^ 2)
    res = isf.grid_res
    velocity = isf.velocity
    map_idx!(isf.f, ()) do xyz, f, _
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
    return
end

function poisson_solve(isf)
    fc = fft(isf.f)
    fc .= (*).(fc, isf.fac)
    ifft!(fc)
    fc
end

function pressure_project!(isf, psi)
    velocity_one_form!(isf, psi)
    div!(isf)
    q = poisson_solve(isf)
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

function schroedinger_flow!(isf, psi)
    # extract single psi values into tmp arrays stored in obj
    psi1 = map(first, psi)
    psi2 = map(last, psi)
    psi1 = fftshift(fft(psi1)); psi2 = fftshift(fft(psi2));
    psi1 .= (*).(psi1, isf.mask)
    psi2 .= (*).(psi2, isf.mask)
    psi1 = ifft!(fftshift(psi1)); psi2 = ifft!(fftshift(psi2));
    psi .= twotuple.(psi1, psi2) # convert back
    psi
end

# simple statically allocated particle type
type Particles{FloatType, IntType}
    xyz::Vector{Point{3, FloatType}}
    active::Vector{IntType}
end
Base.length(p::Particles) = length(p.active)

function staggered_advect!(particle, isf)
    velocity, dt, d, res = isf.velocity, isf.dt, isf.d, isf.grid_res
    gridsize = isf.physical_size
    map!(particle, particle) do p

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
    p   = mod.(Vec(point), gs)
    i   = Vec{3, Int}(floor.(p ./ d)) + 1
    ip  = mod.(i, res) + 1

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
    inserted = 0
    for i = 1:length(p.xyz) # reuse freed up indices
        if !(i in p.active)
            inserted += 1
            p.xyz[i] = xyz[inserted]
            push!(p.active, i)
        end
        inserted == length(xyz) && return
    end
    sort!(p.active)
    length(xyz) == inserted && return
    # We have more particles than `gaps`, meaning all gaps should be filled
    # and length(p.active) == last(p.active)
    @show inserted length(xyz)
    @assert isempty(p.active) || (length(p.active) == last(p.active)) # lets just make sure of this
    rest = length(xyz) - inserted
    newlen = length(p.active) + rest
    if newlen > length(p.xyz)
        # we don't grow particles for now. Statically allocated, Remember?
        warn("Not all particles are inserted, limit reached")
        newlen = length(p.xyz)
    end
    range = (length(p.active)+1):newlen
    p.xyz[range] = view(xyz, (inserted+1):length(xyz))
    append!(p.active, range)
    return
end


end
