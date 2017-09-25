using SchroedingersSmoke, CLArrays
using SchroedingersSmoke.ParallelPort
using Colors, GPUArrays

vol_size = (4,2,2)# box size
dims = (64, 32, 32) .* 2 # volume resolution
hbar = 0.1f0      # Planck constant
dt = 1f0/48f0     # time step

jet_velocity = (1f0, 0f0, 0f0)
nozzle_cen = Float32.((2-1.7, 1-0.034, 1+0.066))
nozzle_len = 0.5f0
nozzle_rad = 0.5f0
n_particles = 500   # number of particles

ArrayType = CLArray

isf2 = ISF{ArrayType, UInt32, Float32}(vol_size, dims, hbar, dt);

tuple_dot(a, b) = sum(a .+ b)

# function returning true at nozzle position
function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end

function restrict_kernel(psi, isjet, kvec, pos, omgterm)
    if isjet
        amp = abs.(psi)
        phase = tuple_dot(kvec, pos) - omgterm
        @fastmath amp .* exp(Complex64(0f0, 1f0) .* phase)
    else
        psi
    end
end

function restrict_velocity!(isf, psi, kvec, isjetarr, omgterm = 0f0)
    psi .= restrict_kernel.(psi, isjetarr, (kvec,), isf.positions, omgterm)
end

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


particles = ArrayType(map(x-> (0f0, 0f0, 0f0), 1:100_000))


newp = map(1:n_particles) do _
    rt = rand()*2*pi
    Float32.((
        nozzle_cen[1] - 0.1,
        nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
        nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
    ))
end;
particles[1:n_particles] = newp

using GLVisualize;
w = glscreen(color = RGBA(0f0, 0f0, 0f0, 0f0));
@async renderloop(w)
particle_vis = visualize(
    (GeometryTypes.Circle(GeometryTypes.Point2f0(0), 0.002f0), GeometryTypes.Point3f0.(Array(particles))),
    boundingbox = nothing, # don't waste time on bb computation
    color = collect(linspace(RGBA{Float32}(1,1,1,0.4), RGBA{Float32}(((243, 123, 173)./255)...,0.4), length(particles))),
    billboard = true
).children[]
_view(particle_vis, camera = :perspective)

# velocity_vis = glplot(GeometryTypes.Vec3f0.(isf.velocity)).children[]

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

        newp = map(1:n_particles) do _
            rt = rand()*2*pi
            offset = ((rand(Float32), rand(Float32), rand(Float32)) .* 0.01f0) .- 0.005f0
            Float32.((
                nozzle_cen[1] - 0.1,
                nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
                nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
            )) .+ offset
        end
        start = mod((iter - 1) * n_particles + 1, length(particles))
        stop = start + n_particles - 1
        particle[start:stop] = newp
        broadcast!(particle, particle, (isf.physical_size,)) do p, ps
            if p[1] > 0.1f0 && p[1] < ps[1] ||
                    p[2] > 0.1f0 && p[2] < ps[2] ||
                    p[3] > 0.1f0 && p[3] < ps[3]
                p
            else
                (0f0, 0f0, 0f0)
            end
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
    200, isf2, psi, kvec, omega, n_particles, isjetarr,
    nozzle_rad, nozzle_cen, particles, particle_vis
)
