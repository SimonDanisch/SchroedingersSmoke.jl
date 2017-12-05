using SchroedingersSmoke, CLArrays
using Colors, GPUArrays

vol_size = (4,2,2)# box size
dims = (64, 32, 32) .* 2 # volume resolution
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


particles = ArrayType(map(x-> (0f0, 0f0, 0f0), 1:(10^6) * 3))

add_particles!(particles, 1:n_particles, nozzle_cen, nozzle_rad)


using GLVisualize;
w = glscreen(color = RGBA(0f0, 0f0, 0f0, 0f0));
@async renderloop(w)
particle_vis = visualize(
    (GeometryTypes.Circle(GeometryTypes.Point2f0(0), 0.002f0), GeometryTypes.Point3f0.(Array(particles))),
    boundingbox = nothing, # don't waste time on bb computation
    color = fill(RGBA{Float32}(0, 0, 0, 0.09), length(particles)),
    billboard = true
).children[]
_view(particle_vis, camera = :perspective)
particle_vis[:color][1:n_particles] = map(1:n_particles) do i
    xx = (i / n_particles) * 2pi
    RGBA{Float32}((sin(xx) + 1) / 2, (cos(xx) + 1.0) / 2.0, 0.0, 0.1)
end
function simloop(
        N, isf, psi, kvec, omega, n_particles, isjetarr,
        nozzle_rad, nozzle_cen, particle, particle_vis
    )
    dt = isf.dt; d = isf.d
    for iter = 1:N
        isopen(w) || break
        t = iter * dt
        # incompressible Schroedinger flow
        schroedinger_flow!(isf, psi)
        psi .= normalize_psi.(psi)
        pressure_project!(isf, psi)

        start = mod((iter - 1) * n_particles + 1, length(particle))
        stop = start + n_particles - 1
        add_particles!(particle, start:stop, nozzle_cen, nozzle_rad)
        particle_vis[:color][start:stop] = map(1:n_particles) do i
            xx = (i / n_particles) * 2pi
            RGBA{Float32}((sin(xx) + 1) / 2, (cos(xx) + 1.0) / 2.0, iter / N, 0.1)
        end
        # constrain velocity
        restrict_velocity!(isf, psi, kvec, isjetarr, omega*t)
        pressure_project!(isf, psi)
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!((x, y)-> x .* y, isf.velocity, isf.velocity, (dinv,))
        staggered_advect!(particle, isf)

        #TODO only update active particles
        GLAbstraction.set_arg!(particle_vis, :position, reinterpret(GeometryTypes.Point3f0, Array(particle)))

        yield()
    end
end

@time simloop(
    300, isf2, psi, kvec, omega, n_particles, isjetarr,
    nozzle_rad, nozzle_cen, particles, particle_vis
)
