Base.FFTW.set_num_threads(8)

using SchroedingersSmoke
import SchroedingersSmoke.ParallelPort
using StaticArrays, Colors, GPUArrays
import GeometryTypes
import ParallelPort: ISF, normalize_psi, pressure_project!
import ParallelPort: velocity_one_form!, schroedinger_flow!
import ParallelPort: Particles, staggered_advect!
import ParallelPort: Vec, Vec3f0, Point, Point3f0
using BenchmarkTools

vol_size = (4,2,2)# box size
dims = (64,32,32) # volume resolution
hbar = 0.1f0      # Planck constant
dt = 1f0/48f0     # time step

jet_velocity = Vec3f0(1, 0, 0)
nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
nozzle_len = 0.5f0
nozzle_rad = 0.5f0
n_particles = 300   # number of particles

isf2 = ISF{Int, Float32}(vol_size, dims, hbar, dt);

# function returning true at nozzle position
function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end


function restrict_velocity!(isf, psi, kvec, isjetarr, omgterm = 0f0)
    broadcast!(psi, psi, isjetarr, (kvec,), isf.positions) do psi, isjet, kvec, pos
        if isjet
            amp = abs.(psi)
            phase = dot(kvec, pos) - omgterm
            amp .* exp(1f0*im .* phase)
        else
            psi
        end
    end
end

# initialize psi
psi = CLArray([(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]);

psi .= normalize_psi.(psi);

kvec = jet_velocity ./ hbar;
omega = sum(jet_velocity.^2f0) / (2f0*hbar);
isjetarr = isjet.(isf2.positions, (nozzle_cen,), (nozzle_len,), (nozzle_rad,))
# constrain velocity
for iter = 1:10
    restrict_velocity!(isf2, psi, kvec, isjetarr, 0f0)
    pressure_project!(isf2, psi)
end

# psi2 ≈ last.(psi)#

## SET PARTICLES
# particle = Particles{Float32, Int}(
#     JLArray(zeros(Point3f0, 100_000)), JLArray(Int[])
# )
particle = Particles{Float32, Int}(
    (zeros(Point3f0, 100_000)), (Int[])
);

function in_grid(i, particle = particle, isf = isf2)
    p = particle.xyz[i]
    for i=1:3
        p[i] > 0.1 && p[i] < isf.physical_size[i] || return false
    end
    true
end
newp = map(1:n_particles) do _
    rt = rand()*2*pi
    Point3f0(
        nozzle_cen[1] - 0.1,
        nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
        nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
    )
end;
append!(particle, newp);
using GLVisualize; w = glscreen(); @async renderloop(w)

particle_vis = visualize(
    (GeometryTypes.Circle(GeometryTypes.Point2f0(0), 0.006f0), GeometryTypes.Point3f0.(particle.xyz)),
    boundingbox = nothing, # don't waste time on bb computation
    color = RGBA{Float32}(0,0,0,0.4),
    billboard = true, indices = particle.active
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
            offset = (rand(Point3f0) .* 0.01f0) - 0.005f0
            Point3f0(
                nozzle_cen[1] - 0.1,
                nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
                nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
            ) .+ offset
        end
        append!(particle, newp)
        filter!(in_grid, particle.active)


        # constrain velocity
        restrict_velocity!(isf, psi, kvec, isjetarr, omega*t)
        pressure_project!(isf, psi)
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!((x, y)-> x .* y, isf.velocity, isf.velocity, (dinv,))
        staggered_advect!(view(particle.xyz, particle.active), isf)

        #TODO only update active particles
        GLAbstraction.set_arg!(particle_vis, :position, reinterpret(GeometryTypes.Point3f0, particle.xyz))
        GLAbstraction.set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active))
        yield()
    end
end


@time simloop(
    200, isf2, psi, kvec, omega, n_particles, isjetarr,
    nozzle_rad, nozzle_cen, particle, particle_vis
)
