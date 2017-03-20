Base.FFTW.set_num_threads(8)

using SchroedingersSmoke
import SchroedingersSmoke.ParallelPort
using StaticArrays, Colors
import GeometryTypes

import ParallelPort: ISF, normalize_psi, pressure_project!
import ParallelPort: velocity_one_form!, schroedinger_flow!
import ParallelPort: Particles, staggered_advect!, map_idx!
import ParallelPort: Vec, Vec3f0, Point, Point3f0, JLArray
using Sugar, BenchmarkTools

vol_size = (4,2,2)# box size
dims = (64,32,32) # volume resolution
hbar = 0.1f0      # Planck constant
dt = 1f0/48f0     # time step

jet_velocity = Vec3f0(1, 0, 0)
nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
nozzle_len = 0.5f0
nozzle_rad = 0.5f0
n_particles = 50   # number of particles

isf = ISF{Int, Float32}(vol_size, dims, hbar, dt);

# function returning true at nozzle position
function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end
function restrict_velocity(pos, psi, args)
    omgterm, kvec, nozzle_cen, nozzle_len, nozzle_rad = args
    if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
        phase = sum(kvec .* pos) - omgterm
        (abs(psi[1]) * exp(1im * phase), abs(psi[2]) * exp(1im * phase))
    end
    psi
end
function restrict_velocity!(isf, psi, kvec, nozzle_cen, nozzle_len, nozzle_rad, omgterm = 1f0)
    args = (omgterm, kvec, nozzle_cen, nozzle_len, nozzle_rad)
    psi .= restrict_velocity.(
        isf.positions,
        psi,
        Scalar(args)
    )
end

# initialize psi
psi = JLArray([(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]);
normalize_psi.(psi);

kvec = jet_velocity ./ hbar;
omega = sum(jet_velocity.^2f0) / (2f0*hbar);

# constrain velocity
for iter = 1:10
    restrict_velocity!(isf, psi, kvec, nozzle_cen, nozzle_len, nozzle_rad)
    pressure_project!(isf, psi)
end


## SET PARTICLES
# particle = Particles{Float32, Int}(
#     JLArray(zeros(Point3f0, 100_000)), JLArray(Int[])
# )
particle = Particles{Float32, Int}(
    (zeros(Point3f0, 100_000)), (Int[])
);

function in_grid(i, particle = particle, isf = isf)
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
using GLPlot; w = GLPlot.init()
particle_vis = glplot(
    (GeometryTypes.Circle(GeometryTypes.Point3f0(0), 0.005f0), GeometryTypes.Point3f0.(particle.xyz)),
    boundingbox = nothing, # don't waste time on bb computation
    color = RGBA{Float32}(0,0,0,0.4),
    billboard = true, indices = particle.active
).children[]

velocity_vis = glplot(GeometryTypes.Vec3f0.(isf.velocity)).children[]
mean(map(first, isf.velocity.buffer))
function simloop(
        N, isf, psi, kvec, omega, n_particles,
        nozzle_rad, nozzle_cen, particle, particle_vis
    )
    dt = isf.dt; d = isf.d
    for iter = 1:N
        t = iter * dt
        # incompressible Schroedinger flow
        schroedinger_flow!(isf, psi)
        normalize_psi.(psi)
        pressure_project!(isf, psi)

        newp = map(1:n_particles) do _
            rt = rand()*2*pi
            Point3f0(
                nozzle_cen[1] - 0.1,
                nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
                nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
            )
        end
        append!(particle, newp)
        filter!(in_grid, particle.active.buffer)

        # constrain velocity
        restrict_velocity!(isf, psi, kvec, nozzle_cen, nozzle_len, nozzle_rad, omega*t)
        pressure_project!(isf, psi)
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!(x-> x .* dinv, isf.velocity, isf.velocity)
        staggered_advect!(view(particle.xyz.buffer, particle.active.buffer), isf)

        #TODO only update active particles
        GLAbstraction.set_arg!(particle_vis, :position, reinterpret(GeometryTypes.Point3f0, particle.xyz.buffer))
        GLAbstraction.set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active.buffer))
        yield()
    end
end
@time simloop(
    200, isf, psi, kvec, omega, n_particles,
    nozzle_rad, nozzle_cen, particle, particle_vis
)


#_view(velocity_vis)

function simloop(
        N, isf, psi, kvec, omega, n_particles,
        nozzle_rad, nozzle_cen, particle, particle_vis, velocity_vis
    )
    dt = isf.dt; d = isf.d
    for iter = 1:N
        t = iter * dt
        # incompressible Schroedinger flow
        schroedinger_flow!(isf, psi)
        normalize_psi.(psi)
        pressure_project!(isf, psi)

        # constrain velocity
        restrict_velocity!(isf, psi, kvec, nozzle_cen, nozzle_len, nozzle_rad, omega*t)
        pressure_project!(isf, psi)

        #set_arg!(velocity_vis, :rotation, vec(isf.velocity))
        # particle birth


        # advect and show particles
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        isf.velocity .= (.*).(isf.velocity, Scalar(inv.(d)))

        staggered_advect!(view(particle.xyz, particle.active), isf)

        #TODO only update active particles
        set_arg!(particle_vis, :position, particle.xyz)
        set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active))
        yield()
    end
end
