const Vec{N, T} = NTuple{N, T}
const Point{N, T} = NTuple{N, T}

struct ISF{
        IntType, FloatType,
        # can't define the Types as e.g. ArrayType{FloatType, 3}, so we need them as type parameter and calculate the types in the constructor
        CArray, FArray, VIArray, VFArray, # Complex array, Float array, Integer Vec array, Float Vec array
        FFTP, IFFTP # types of the fft plan
    }

    grid_res::Vec{3, IntType}
    physical_size::Vec{3, IntType}
    d::Vec{3, FloatType}
    velocity::VFArray

    hbar::FloatType             # reduced Planck constant
    dt::FloatType               # time step
    mask::CArray # Fourier coefficient for solving Schroedinger eq
    # precalculated values
    # we have a couple of indices precalculated:
    i::VIArray # type is to complex to write down
    i_shifted::VIArray # indices shifted by one
    i_shifted2::VIArray # indices shifted by -2

    # other precalculated values
    positions::VFArray # position on the grid
    fac::CArray
    # temporary allocations for intermediates
    f::FArray
    fc::CArray
    psi1::CArray
    psi2::CArray
    fftplan::FFTP
    ifftplan::IFFTP
end

function Base.zero(::Type{NTuple{N, T}}) where {T, N}
    ntuple(x-> zero(T), Val{N})
end

function (::Type{ISF{ArrayType, IntType, FloatType}})(physical_size, dims, hbar, dt) where {ArrayType, IntType, FloatType}
    VF = Vec{3, FloatType}
    VI = Vec{3, IntType}
    grid_res = VI(dims)
    physical_size = VI(physical_size)
    d = VF(physical_size ./ grid_res)
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
    velocity = zeros(VF, dims)

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
    ISF{
            IntType, FloatType,
            ArrayType{Complex{FloatType}, 3}, ArrayType{FloatType, 3},
            ArrayType{VI, 3}, ArrayType{VF, 3},
            typeof(p1), typeof(p2)
        }(
        grid_res, physical_size, d,
        velocity,
        hbar, dt,
        mask,
        i, i_shifted, i_shifted2,
        positions,
        fac_gpu, f,
        copy(fac_gpu), # temporaries
        copy(fac_gpu),
        copy(fac_gpu),
        p1, p2
    )
end

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
    norm = @fastmath hypot(abs(psi[1]), abs(psi[2]))
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


function in_grid(p, physical_size)
    p[1] > 0.0f0 && p[1] < physical_size[1] || return false
    p[2] > 0.0f0 && p[2] < physical_size[2] || return false
    p[3] > 0.0f0 && p[3] < physical_size[3] || return false
    true
end

function staggered_advect(p, velocity, args)
    dt = args[1]; d = args[2]; gridsize = args[3]; res = args[4]
    if in_grid(p, gridsize)

        k1 = staggered_velocity(velocity, p, d, gridsize, res)

        k2 = p .+ k1 .* dt .* 0.5f0
        k2 = staggered_velocity(velocity, k2, d, gridsize, res)

        k3 = p .+ k2 .* dt .* 0.5f0
        k3 = staggered_velocity(velocity, k3, d, gridsize, res)

        k4 = p .+ k3 .* dt
        k4 = staggered_velocity(velocity, k4, d, gridsize, res)

        p .+ dt./6f0 .* (k1 .+ 2f0 .* k2 .+ 2f0 .* k3 .+ k4)
    else
        (-1f0, -1f0, -1f0)
    end
end


@inline function staggered_velocity(velocity, point, d, ph_size, res)
    p   = mod.(point, ph_size)
    i   = UInt32.(floor.(p ./ d)) .+ UInt32(1)
    ip  = mod.(i, res) .+ UInt32(1)

    v0  = velocity[i[1], i[2], i[3]]

    pxp = velocity[ip[1], i[2], i[3]]
    pyp = velocity[i[1], ip[2], i[3]]
    pzp = velocity[i[1], i[2], ip[3]]

    vn = (
        velocity[i[1], ip[2], ip[3]][1],
        velocity[ip[1], i[2], ip[3]][2],
        velocity[ip[1], ip[2], i[3]][3]
    )
    pp  = (pyp[1], pxp[2], pxp[3])
    pp2 = (pzp[1], pzp[2], pyp[3])

    w   = p .- (i .- 1f0) .* d
    w1  = (w[3], w[3], w[2])
    w2  = (w[2], w[1], w[1])

    return (
        (1f0 .- w1) .* ((1f0 .- w2) .* v0 .+ w2 .* pp) .+
        w1 .* ((1f0 .- w2) .* pp2 .+ w2 .* vn)
    )
end

function staggered_advect!(particle, isf)
    velocity, dt, d, res = isf.velocity, isf.dt, isf.d, isf.grid_res
    gridsize = isf.physical_size
    particle .= staggered_advect.(
        particle,
        (velocity,),
        ((
            isf.dt,
            isf.d,
            gridsize,
            isf.grid_res
        ),)
    )
end


export ISF, normalize_psi, pressure_project!
export velocity_one_form!, schroedinger_flow!
export Particles, staggered_advect!
