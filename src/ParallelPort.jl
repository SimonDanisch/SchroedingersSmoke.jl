module ParallelPort

using StaticArrays, CLArrays
const Vec{N, T} = NTuple{N, T}
const Point{N, T} = NTuple{N, T}
const Vec3f0 = NTuple{3, Float32}
const Point3f0 = NTuple{3, Float32}

const ArrayType = CLArray

# lots of parameters. To lazy to write out types,
# still don't want to waste performance
type ISF{IntType, FloatType, FFTP, IFFTP}

    grid_res::Vec{3, IntType}
    physical_size::Vec{3, IntType}
    d::Vec{3, FloatType}
    velocity::ArrayType{Vec{3, FloatType}, 3}

    hbar::FloatType             # reduced Planck constant
    dt::FloatType               # time step
    mask::ArrayType{Complex{FloatType}, 3} # Fourier coefficient for solving Schroedinger eq
    # precalculated values
    # we have a couple of indices precalculated:
    i::ArrayType{NTuple{3, IntType}, 3} # type is to complex to write down
    i_shifted::ArrayType{Vec{3, IntType}, 3} # indices shifted by one
    i_shifted2::ArrayType{Vec{3, IntType}, 3} # indices shifted by -2

    # other precalculated values
    positions::ArrayType{Point{3, FloatType}, 3} # position on the grid
    fac::ArrayType{Complex{FloatType}, 3}
    # temporary allocations for intermediates
    f::ArrayType{FloatType, 3}
    fc::ArrayType{Complex{FloatType}, 3}
    psi1::ArrayType{Complex{FloatType}, 3}
    psi2::ArrayType{Complex{FloatType}, 3}
    fftplan::FFTP
    ifftplan::IFFTP
end

function Base.zero(::Type{NTuple{N, T}}) where {T, N}
    ntuple(x-> zero(T), Val{N})
end

function (::Type{ISF{IntType, FloatType}})(physical_size, dims, hbar, dt) where {IntType, FloatType}
    VT = NTuple{3, FloatType}
    grid_res = Vec(dims)
    physical_size = Vec(physical_size)
    d = VT(physical_size ./ grid_res)
    fac = zeros(Complex{FloatType}, dims)
    mask = zeros(Complex{FloatType}, dims)
    f = -4 * pi^2 * hbar
    for x = 1:dims[1], y = 1:dims[2], z = 1:dims[3]
        xyz = (x, y, z)
        # fac
        s = sin.(pi .* (xyz .- 1) ./ grid_res) ./ d
        denom = sum(s .^ 2)
        fac[x,y,z] = -0.25f0 ./ denom

        # schroedingers mask
        k = (xyz .- 1 .- grid_res ./ 2) ./ physical_size
        lambda = f .* sum(k .^ 2)
        mask[x,y,z] = exp(1.0im * lambda * dt / 2)
    end
    mask = fftshift(mask)
    fac[1,1,1] = 0

    f = zeros(FloatType, dims)
    velocity = zeros(VT, dims)

    i = collect(((x, y, z) for x=1:dims[1], y=1:dims[2], z=1:dims[3]))
    i_shifted = map(i) do i
        mod.(Vec(i), Vec(dims)) .+ 1
    end
    i_shifted2 = map(i) do i
         mod.(i .- 2, dims) .+ 1
    end
    positions = map(i) do i
        (i .- 1) .* d
    end
    fac_gpu = ArrayType(fac)
    p1, p2 = plan_fft!(fac_gpu), plan_ifft!(fac_gpu)
    ISF{IntType, FloatType, typeof(p1), typeof(p2)}(
        grid_res,
        physical_size,
        d,
        velocity,
        hbar,
        dt,
        mask,
        i, i_shifted, i_shifted2,
        positions,
        fac,
        f,
        complex.(f),
        complex.(f),
        complex.(f),
        p1, p2
    )
end

# helper function to have a map but with random index manipulations

@inline function velocity_one_form(i, i2, psi, hbar)
    @inbounds begin
        psi12  = psi[i[1],  i[2],  i[3]]
        psix12 = psi[i2[1], i[2],  i[3]]
        psiy12 = psi[i[1],  i2[2] ,i[3]]
        psiz12 = psi[i[1],  i[2],  i2[3]]
        psi1n = (psix12[1], psiy12[1], psiz12[1])
        psi2n = (psix12[2], psiy12[2], psiz12[2])
        return angle.(
            conj(psi12[1]) .* psi1n .+
            conj(psi12[2]) .* psi2n
        ) .* hbar
    end
end
function velocity_one_form!(isf, psi, hbar = 1.0f0)
    isf.velocity .= velocity_one_form.(
        isf.i,
        isf.i_shifted,
        (psi,),
        hbar
    )
end



@fastmath function gauge_transform(psi, q)
    x = Complex64(0f0, 1f0) * (-q)
    eiq = exp(x)
    (psi[1] * eiq, psi[2] * eiq)
end

function div(i, i2, velocity, res, ds)
    x = i[1]; y = i[2]; z = i[3]
    ix = i2[1]; iy = i2[2]; iz = i2[3]
    v1 = velocity[x, y,  z]
    v2 = (
        velocity[ix, y,  z][1],
        velocity[x,  iy, z][2],
        velocity[x,  y, iz][3]
    )
    t1 = v1 .- v2
    t2 = t1 .* ds
    sum(t2)
end
function div!(isf)
    isf.f .= div.(
        isf.i,
        isf.i_shifted2,
        (isf.velocity,),
        (isf.grid_res,),
        (inv.(isf.d .^ 2f0),)
    )
end


function pressure_project!(isf, psi)
    velocity_one_form!(isf, psi)
    div!(isf)
    q = poisson_solve(isf)
    psi .= gauge_transform.(psi, q)
end

function normalize_psi(psi)
    norm = hypot(abs(psi[1]), abs(psi[2]))
    (psi[1] / norm, psi[2] / norm)
end
@inline twotuple(a, b) = (a, b)

function poisson_solve(isf)
    isf.fc .= complex.(isf.f)
    isf.fftplan * isf.fc
    isf.fc .= (*).(isf.fc, isf.fac)
    isf.ifftplan * isf.fc
    isf.fc
end

function schroedinger_flow!(isf, psi)
    # extract single psi values into tmp arrays stored in obj
    psi1 = isf.psi1; psi2 = isf.psi2
    psi1 .= getindex.(psi, 1)
    psi2 .= getindex.(psi, 2)
    isf.fftplan * psi1; isf.fftplan * psi2;
    psi1 .= (*).(psi1, isf.mask)
    psi2 .= (*).(psi2, isf.mask)
    isf.ifftplan * psi1; isf.ifftplan * psi2;
    psi .= twotuple.(psi1, psi2) # convert back
    psi
end

# simple statically allocated particle type
type Particles{FloatType, IntType}
    xyz::ArrayType{Point{3, FloatType}, 1}
    active::ArrayType{IntType, 1}
end
Base.length(p::Particles) = length(p.active)


function staggered_advect(p, args)
    velocity, dt, d, gridsize, res = args
    k1 = staggered_velocity(velocity, p, d, gridsize, res)

    k2 = p + k1 .* dt * 0.5f0
    k2 = staggered_velocity(velocity, k2, d, gridsize, res)

    k3 = p + k2 .* dt * 0.5f0
    k3 = staggered_velocity(velocity, k3, d, gridsize, res)

    k4 = p + k3 .* dt
    k4 = staggered_velocity(velocity, k4, d, gridsize, res)

    p .+ dt/6f0 .* (k1 .+ 2f0*k2 .+ 2f0*k3 .+ k4)
end
function staggered_advect!(particle, isf)
    velocity, dt, d, res = isf.velocity, isf.dt, isf.d, isf.grid_res
    gridsize = isf.physical_size
    particle .= staggered_advect.(
        particle,
        ((
            velocity,
            isf.dt,
            isf.d,
            isf.physical_size,
            isf.grid_res
        ),)
    )
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
