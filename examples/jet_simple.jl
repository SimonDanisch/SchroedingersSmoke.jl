Base.FFTW.set_num_threads(8)
BLAS.set_num_threads(8)

using SchroedingersSmoke

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



# Set nozzle
isJet = land(
    abs(isf.t.px - nozzle_cen[1]).<=nozzle_len/2,
    (isf.t.py - nozzle_cen[2]).^2+(isf.t.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
)

# initialize psi
psi1f = ones(size(isf.t.px))
psi2f = psi1f * 0.01
psi1f, psi2f = Normalize(psi1f, psi2f)

# constrain velocity
kvec = jet_velocity / isf.hbar
omega = sum(jet_velocity .^ 2) / (2*isf.hbar);
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;
# convert to complex
psi1 = (1.+0.0im)*psi1f
psi2 = (1.+0.0im)*psi2f
for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet] .* exp(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet] .* exp(1.0im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end
max_particles = 20_000
max_history = 19

function removethedead!(particles, vol_size)
    map!(particles.active) do i
        !isinside(
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
    0, active_particles
)


function loop(
        N, isf, psi1, psi2, dt, kvec, omega, isJet, n_particles,
        nozzle_rad, nozzle_cen, update,
        particle, gpu_position_x, gpu_position_y, gpu_position_z,
    )
    for iter = 1:N
        t = iter * dt
        # incompressible Schroedinger flow
        psi1, psi2 = SchroedingerFlow(isf, psi1, psi2)
        psi1, psi2 = Normalize(psi1, psi2)
        psi1, psi2 = PressureProject(isf, psi1, psi2)

        # constrain velocity
        phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz - omega*t
        amp1 = abs(psi1)
        amp2 = abs(psi2)
        psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
        psi1, psi2  = PressureProject(isf, psi1, psi2)

        # particle birth
        rt   = rand(Float32, n_particles)*2*pi
        newx = nozzle_cen[1] * ones(Float32, size(rt))
        newy = nozzle_cen[2] + 0.9*nozzle_rad*cos(rt)
        newz = nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)

        append!(particle, newx, newy, newz)
        removethedead!(particles, vol_size)
        # advect and show particles
        vx,vy,vz = VelocityOneForm(isf, psi1, psi2, isf.hbar)
        vx,vy,vz = StaggeredSharp(isf.t, vx,vy,vz)
        StaggeredAdvect(particle, isf.t, vx,vy,vz, isf.dt)

    end
end
@time loop(
    20, isf, psi1, psi2, dt, kvec, omega, isJet, n_particles,
    nozzle_rad, nozzle_cen,
)
# Profile.print()
#
# using ProfileView
# ProfileView.view()
