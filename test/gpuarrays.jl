using GPUArrays, StaticArrays
CLBackend.init()
import Base: RefValue
const Vec = SVector
immutable ISF4{
        IdxT <: AbstractArray{T, 3} where T <: Vec{3, IT} where IT <: Integer,
        VeloT <: AbstractArray{T2, 3} where T2 <: Vec{3, FT} where FT <: AbstractFloat
    }
    i::IdxT
    i_shifted::IdxT
    velocity::VeloT
end
dims = (100, 100, 100)
IntType = UInt32
const IVec2 = Vec{3, IntType}
i = [IVec2(IntType(x), IntType(y), IntType(z)) for x=1:dims[1], y=1:dims[2], z=1:dims[3]]
i_shifted = map(i) do i
    mod.(IVec2(i), IVec2(dims)) + IntType(1)
end
isf = ISF4(GPUArray(i), GPUArray(i_shifted), GPUArray(zeros(Vec{3, Float32}, dims)))
psi = GPUArray(collect(zip(rand(Complex64, dims), rand(Complex64, dims))))
using GPUArrays: gpu_sub2ind


@inline function velocity_one_form(i, i2, psi, hbar)
    s = (Cuint(100), Cuint(100), Cuint(100))
    psi12  = psi[gpu_sub2ind(s, (i[1],  i[2],  i[3]))]
    psix12 = psi[gpu_sub2ind(s, (i2[1], i[2],  i[3]))]
    psiy12 = psi[gpu_sub2ind(s, (i[1],  i2[2] ,i[3]))]
    psiz12 = psi[gpu_sub2ind(s, (i[1],  i[2],  i2[3]))]
    psi1n  = Vec(psix12[1], psiy12[1], psiz12[1])
    psi2n  = Vec(psix12[2], psiy12[2], psiz12[2])
    angle.(
        conj(psi12[1]) * psi1n .+
        conj(psi12[2]) * psi2n
    ) * hbar
end

function velocity_one_form!(isf, psi, hbar = 1.0f0)
    isf.velocity .= velocity_one_form.(
        isf.i,
        isf.i_shifted,
        RefValue(psi),
        hbar
    )
end
velocity_one_form!(isf, psi)

i = [IVec2(IntType(x), IntType(y), IntType(z)) for x=1:dims[1], y=1:dims[2], z=1:dims[3]]
i_shifted = map(i) do i
    mod.(IVec2(i), IVec2(dims)) + IntType(1)
end
isf = ISF4((i), (i_shifted), (zeros(Vec{3, Float32}, dims)))
psi = (collect(zip(rand(Complex64, dims), rand(Complex64, dims))))
velocity_one_form!(isf, psi)

x = GPUArray(NTuple{3, Cuint}[(0, 0, 0)])
y = GPUArray(reinterpret(Float32, [(NTuple{3, Cuint}((1, 2, 3))..., 1f0, 2f0, 3f0)]))

function kernel(state, x, y)
    x[1] = Transpiler.cli.CLArray{NTuple{3, Cuint}, 1}(y)[1]
    return
end
println(Array(x))
gpu_call(kernel, x, (x, y))
Array(x)
