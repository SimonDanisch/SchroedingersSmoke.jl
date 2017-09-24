using SchroedingersSmoke, CLArrays
import SchroedingersSmoke.ParallelPort
using StaticArrays, Colors

import ParallelPort: ISF, normalize_psi, pressure_project!
import ParallelPort: velocity_one_form!, schroedinger_flow!
import ParallelPort: Particles, staggered_advect!

# function returning true at nozzle position
function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end

Base.dot(a::NTuple{N, T}, b::NTuple{N, T}) where {N, T} = sum(a .* b)
Complex64(0f0, 1f0) ==  1f0*im

function restrict_kernel(psi, isjet, kvec, pos, omgterm)
    if isjet
        amp = abs.(psi)
        phase = dot(kvec, pos) - omgterm
        @fastmath amp .* exp(Complex64(0f0, 1f0) .* phase)
    else
        psi
    end
end

function restrict_velocity!(isf, psi, kvec, isjetarr, omgterm = 0f0)
    psi .= restrict_kernel.(psi, isjetarr, (kvec,), isf.positions, omgterm)
end


vol_size = (4,2,2)# box size
dims = (64,32,32) # volume resolution
hbar = 0.1f0      # Planck constant
dt = 1f0/48f0     # time step

jet_velocity = Float32.((1, 0, 0))
nozzle_cen = Float32.((2-1.7, 1-0.034, 1+0.066))
nozzle_len = 0.5f0
nozzle_rad = 0.5f0
n_particles = 50   # number of particles

isf = ISF{UInt32, Float32}(vol_size, dims, hbar, dt);


# initialize psi
psi = CLArray([(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]);
psi .= normalize_psi.(psi);
kvec = jet_velocity ./ hbar;
omega = sum(jet_velocity.^2f0) / (2f0*hbar);
isjetarr = isjet.(isf.positions, (nozzle_cen,), (nozzle_len,), (nozzle_rad,))

for iter = 1:10
    restrict_velocity!(isf, psi, kvec, isjetarr, 0f0)
    pressure_project!(isf, psi)
end

function simloop(
        N, isf, psi, kvec, omega,
        isjetarr
    )
    dt = isf.dt; d = isf.d
    for iter = 1:N
        t = Float32(iter * dt)
        # incompressible Schroedinger flow
        schroedinger_flow!(isf, psi)
        psi .= normalize_psi.(psi)
        pressure_project!(isf, psi)

        # constrain velocity
        restrict_velocity!(isf, psi, kvec, isjetarr, omega*t)
        pressure_project!(isf, psi)
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!((a, b)-> a .* b, isf.velocity, isf.velocity, (dinv,))
        #staggered_advect!(view(particle.xyz, particle.active), isf)
    end
end

@time begin
    simloop(200, isf, psi, kvec, omega, isjetarr)
    GPUArrays.synchronize(psi)
end






particle = CLArray(map(x-> (rand(Float32), rand(Float32), rand(Float32)) .* (4f0, 2f0, 2f0), 1:10_000))
@time begin
    staggered_advect!(particle, isf)
    GPUArrays.synchronize(particle)
end

gpu_call(particle, (velocity, isf.dt, isf.d, isf.physical_size, isf.grid_res)) do state, velo, dt, d, ps, gr
    point = (0.5f0, 0.5f0, 0.5f0)
    staggered_velocity(velo, point, d, ps, gr)
    return
end

typename(CL)
particle .= staggered_advect.(
    particle,
    (isf.velocity,),
    ((
        isf.dt,
        isf.d,
        isf.physical_size,
        isf.grid_res
    ),)
)
