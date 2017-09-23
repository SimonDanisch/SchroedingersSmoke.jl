using CLArrays, StaticArrays
import Base: RefValue
import GPUArrays: JLArray


const dims = (100, 100, 100)
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

using CLArrays
Typ = CLArray
i = [(UInt32(x), UInt32(y), UInt32(z)) for x=1:dims[1], y=1:dims[2], z=1:dims[3]]
i_shifted = map(i) do i
    mod.(UInt32.(i), UInt32.(dims)) .+ UInt32(1)
end
velocity = Typ(fill(ntuple(x-> 0f0, Val{3}), dims))
i = Typ(i)
i_shifted = Typ(i_shifted)
psi = Typ(collect(zip(rand(Complex64, dims), rand(Complex64, dims))))
f = similar(psi, Float32)


@fastmath function normalize_psi2(psi)
    norm = hypot(abs(psi[1]), abs(psi[2]))
    (psi[1] / norm, psi[2] / norm)
end

psi .= normalize_psi.(psi)

using BenchmarkTools

a, b, c, d = Array.((velocity, i, i_shifted, psi))

function gauge_transform3(psi, q)
    x = Complex64(0f0, 1f0) * (-q)
    eiq = exp(x)
    (psi[1] * eiq, psi[2] * eiq)
end

q_c = rand(Complex64, dims)

psi_in = collect(zip(q_c, q_c))

qg, psig = Typ(q_c), Typ(psi_in)

@btime begin
    $psi .= normalize_psi2.($psi)
    GPUArrays.synchronize($psi)
end


using CuArrays

function divide(i, i2, velocity, ds)
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
ds = inv.((0.1f0, 0.1f0, 0.1f0) .^ 2f0)
@btime begin
    $f .= $divide.($i, $i_shifted, $((velocity,)), $((ds,)))
    GPUArrays.synchronize($f)
end
@which GPUArrays.synchronize(f)
fc, ic, isc, vc = Array.((f, i,  i_shifted, velocity))
@btime $fc .= $divide.($ic, $isc, $((vc,)), $((ds,)))
