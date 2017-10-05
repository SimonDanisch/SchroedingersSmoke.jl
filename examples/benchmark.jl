using SchroedingersSmoke, CLArrays
using Colors

vol_size = (4,2,2)# box size
dims = (64, 32, 32) .* 4 # volume resolution
hbar = 0.1f0      # Planck constant
dt = 1f0/48f0     # time step

jet_velocity = (1f0, 0f0, 0f0)
nozzle_cen = Float32.((2-1.7, 1-0.034, 1+0.066))
nozzle_len = 0.5f0
nozzle_rad = 0.5f0
n_particles = 1000   # number of particles

ArrayType = CLArray

isf2 = ISF{ArrayType, UInt32, Float32}(vol_size, dims, hbar, dt);

# initialize psi
psi = ArrayType([(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]);

psi .= normalize_psi.(psi);

kvec = jet_velocity ./ hbar;
omega = sum(jet_velocity.^2f0) / (2f0*hbar);
isjetarr = isjet.(isf2.positions, (nozzle_cen,), (nozzle_len,), (nozzle_rad,))
# constrain velocity
for iter = 1:10
    restrict_velocity!(isf2, psi, kvec, isjetarr, 0f0)
    pressure_project!(isf2, psi)
end

# add out of bounds particles
particles = ArrayType(map(x-> (-1f0, -1f0, -1f0), 1:10^6))



add_particles!(particles, 1:n_particles, nozzle_cen, nozzle_rad)


function simloop(
        N, isf, psi, kvec, omega, n_particles, isjetarr,
        nozzle_rad, nozzle_cen, particle
    )
    dt = isf.dt; d = isf.d
    for iter = 1:N
        t = iter * dt
        # incompressible Schroedinger flow
        schroedinger_flow!(isf, psi)
        psi .= normalize_psi.(psi)
        pressure_project!(isf, psi)

        start = mod((iter - 1) * n_particles + 1, length(particles))
        stop = start + n_particles - 1
        add_particles!(particles, start:stop, nozzle_cen, nozzle_rad)

        # constrain velocity
        restrict_velocity!(isf, psi, kvec, isjetarr, omega*t)
        pressure_project!(isf, psi)
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!((x, y)-> x .* y, isf.velocity, isf.velocity, (dinv,))
        staggered_advect!(particle, isf)
    end
    GPUArrays.synchronize(particle)
end


@time simloop(
    50, isf2, psi, kvec, omega, n_particles, isjetarr,
    nozzle_rad, nozzle_cen, particles
)
