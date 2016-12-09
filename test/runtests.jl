using GeometryTypes
function test(velocity, res, psi, hbar)
    @inbounds for z = 1:size(velocity, 3), y = 1:size(velocity, 2)
        @simd for x = 1:size(velocity, 1)
            i2 = mod.(Vec(x,y,z), res) + 1
            psi12  = psi[x,  y,  z]
            psix12 = psi[i2[1], y,  z]
            psiy12 = psi[x,  i2[2] ,z]
            psiz12 = psi[x,  y,  i2[3]]
            psi1n = Vec(psix12[1], psiy12[1], psiz12[1])
            psi2n = Vec(psix12[2], psiy12[2], psiz12[2])
            velocity[x,y,z] = angle.(
                conj(psi12[1]) .* psi1n .+
                conj(psi12[2]) .* psi2n
            ) * hbar
        end
    end
end

dims = (64, 32, 32)
grid_res = Vec(dims)
hbar = 0.1f0; dt = 1f0/48f0;
grid_size = Vec(4, 2, 2)
res = Vec(grid_res)

velocity = zeros(Vec3f0, dims)
psi = [(one(Complex64), one(Complex64) * 0.1f0)
for i=1:dims[1], j=1:dims[2], k=1:dims[3]]

using BenchmarkTools

psi1 = map(first, psi)
psi2 = map(last, psi)
@time test(velocity, res, psi, 1f0)
