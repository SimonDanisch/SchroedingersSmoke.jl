Base.FFTW.set_num_threads(8)
BLAS.set_num_threads(8)
include("../../src/MatlabPort.jl")

using GeometryTypes, MatlabPort, GLAbstraction
import MatlabPort: ISF, Normalize, PressureProject, land, Particles, StaggeredAdvect
import MatlabPort: VelocityOneForm, StaggeredSharp, SchroedingerFlow

function GeometryTypes.isinside(point, lower, upper)
    @inbounds for (p,l,u) in zip(point, lower, upper)
        (p > l && p < u) || return false
    end
    true
end


vol_size = (4,2,2)   # box size
vol_res = (64,32,32) # volume resolution
hbar = 0.1           # Planck constant
dt = 1/48            # time step
tmax = 50            # max time

jet_velocity = [1,0,0]; # jet velocity

nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; # nozzle center
nozzle_len = 0.5;                   # nozzle length
nozzle_rad = 0.5;                   # nozzle radius

n_particles = 10;   # number of particles

isf = MatlabPort.ISF(vol_size, vol_res, hbar, dt)

# Set nozzle
isJet = land(
    abs.(isf.px .- nozzle_cen[1]) .<= nozzle_len./2,
    (isf.py .- nozzle_cen[2]) .^ 2 .+ (isf.pz .- nozzle_cen[3]) .^ 2 .<= nozzle_rad .^ 2
)

# initialize psi
psi1f = ones(size(isf.px))
psi2f = psi1f * 0.01
psi1f, psi2f = Normalize(psi1f, psi2f)

# constrain velocity
kvec = jet_velocity / isf.hbar
omega = sum(jet_velocity .^ 2) ./ (2*isf.hbar);
phase = kvec[1].*isf.px .+ kvec[2].*isf.py .+ kvec[3].*isf.pz;
# convert to complex
psi1 = (1.+0.0im) .* psi1f
psi2 = (1.+0.0im) .* psi2f
for iter = 1:10
    amp1 = abs.(psi1)
    amp2 = abs.(psi2)
    psi1[isJet] = amp1[isJet] .* exp.(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet] .* exp.(1.0im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end


max_particles = 20_000
max_history = 19

function removethedead!(particles, vol_size)
    filter!(particles.active) do i
        isinside(
            (particles.x[i], particles.y[i], particles.z[i]),
            (0,0,0), vol_size
        )
    end
end


## SET PARTICLES
particle = Particles(
    zeros(Float32, max_particles),
    zeros(Float32, max_particles),
    zeros(Float32, max_particles),
    0, Int[]
)


using GLVisualize; w = glscreen(); @async GLWindow.renderloop(w)

particle_vis = visualize(
    (Circle(Point2f0(0), 0.006f0), (particle.x, particle.y, particle.z)),
    boundingbox = nothing, # don't waste time on bb computation
    billboard = true, indices = particle.active
).children[]
_view(particle_vis, camera = :perspective)


function loop(
        N, isf, psi1, psi2, dt, kvec, omega, isJet, n_particles,
        nozzle_rad, nozzle_cen, particle, particle_vis
    )
    gpu_position_x = particle_vis[:position_x]
    gpu_position_y = particle_vis[:position_y]
    gpu_position_z = particle_vis[:position_z]
    for iter = 1:N
        t = iter * dt
        # incompressible Schroedinger flow
        psi1, psi2 = SchroedingerFlow(isf, psi1, psi2)
        psi1, psi2 = Normalize(psi1, psi2)
        psi1, psi2 = PressureProject(isf, psi1, psi2)

        # constrain velocity
        phase = kvec[1] .* isf.px .+ kvec[2] .* isf.py .+ kvec[3] .* isf.pz .- omega .* t
        amp1 = abs.(psi1)
        amp2 = abs.(psi2)
        psi1[isJet] = amp1[isJet] .* exp.(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet] .* exp.(1.0im*phase[isJet])
        psi1, psi2  = PressureProject(isf, psi1, psi2)

        # particle birth
        rt   = rand(Float32, n_particles)*2*pi
        newx = nozzle_cen[1] .* ones(Float32, size(rt))
        newy = nozzle_cen[2] .+ 0.9 .* nozzle_rad .* cos.(rt)
        newz = nozzle_cen[3] .+ 0.9 .* nozzle_rad .* sin.(rt)
        append!(particle, newx, newy, newz)
        removethedead!(particle, vol_size)
        # advect and show particles
        vx,vy,vz = VelocityOneForm(isf, psi1, psi2, isf.hbar)
        vx,vy,vz = StaggeredSharp(isf, vx,vy,vz)
        StaggeredAdvect(particle, isf, vx,vy,vz, isf.dt)

        # update gpu buffer
        r = 1:length(particle)
        gpu_position_x[r] = particle.x[r]
        gpu_position_y[r] = particle.y[r]
        gpu_position_z[r] = particle.z[r]
        set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active))
        yield()
    end
end
@time loop(
    1000, isf, psi1, psi2, dt, kvec, omega, isJet, n_particles,
    nozzle_rad, nozzle_cen, particle, particle_vis
)
