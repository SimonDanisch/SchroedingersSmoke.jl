Base.FFTW.set_num_threads(4)
blas_set_num_threads(4)
using SchroedingersSmoke, GeometryTypes

# example_jet
# An example of incompressible Schroedinger flow producing a jet.
#
## PARAMETERS
const vol_size = (4,2,2);   # box size
const vol_res = (64,32,32); # volume resolution
const hbar = 0.1;            # Planck constant
const dt = 1/48;             # time step
const tmax = 50;             # max time

# another interesting set of parameter:
# const vol_size = (4,2,2);   # box size
# const vol_res = (128,64,64); # volume resolution
# const hbar = 0.02;            # Planck constant
# const dt = 1/48;             # time step
# const tmax = 10;             # max time

const jet_velocity = [1,0,0]; # jet velocity

const nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; # nozzle center
const nozzle_len = 0.5;                   # nozzle length
const nozzle_rad = 0.5;                   # nozzle radius

const n_particles = 50;   # number of particles


## INITIALIZATION


isf = ISF(TorusDEC(vol_size, vol_res), hbar, dt)

# Set nozzle
const isJet = land(
    abs(isf.t.px - nozzle_cen[1]).<=nozzle_len/2,
    (isf.t.py - nozzle_cen[2]).^2+(isf.t.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
)

# initialize psi
psi1f = ones(size(isf.t.px))
psi2f = psi1f*0.01
psi1f, psi2f = Normalize(psi1f, psi2f)

# constrain velocity
kvec = jet_velocity/isf.hbar;
omega = sum(jet_velocity.^2)/(2*isf.hbar);
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;

# convert to complex
psi1 = (1.+0f0*im)*psi1f
psi2 = (1.+0f0*im)*psi2f
for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet].*exp(1f0*im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1f0*im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end


## SET PARTICLES
particle = Particles(Point3f0[])

function GeometryTypes.isinside(point, lower, upper)
    @inbounds for (p,l,u) in zip(point, lower, upper)
        (p > l && p < u) || return false
    end
    true
end
## MAIN ITERATION
itermax = ceil(tmax/dt);
function iterate(particle, isf, psi1, psi2, iter, omega, isJet)
    t = iter*dt
    # incompressible Schroedinger flow
    psi1, psi2 = SchroedingerFlow(isf, psi1,psi2)
    psi1, psi2 = Normalize(psi1,psi2)
    psi1, psi2 = PressureProject(isf, psi1,psi2)

    # constrain velocity
    phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz - omega*t
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet].*exp(1f0*im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1f0*im*phase[isJet])

    psi1, psi2 = PressureProject(isf, psi1, psi2)

    # particle birth
    newp = Point3f0[begin
        rt = rand()*2*pi
        Point3f0(
            nozzle_cen[1],
            nozzle_cen[2] + 0.9*nozzle_rad*cos(rt),
            nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)
        )
    end for i=1:n_particles]

    append!(particle.xyz, newp)
    # advect and show particles
    velocity = VelocityOneForm(isf, psi1, psi2, isf.hbar)
    velocity = StaggeredSharp(isf.t, velocity)
    StaggeredAdvect(particle, isf.t, velocity, isf.dt)

    keep = [isinside(p, Point3f0(0f0), vol_size) for p in particle.xyz]
    Keep(particle, keep)
end

for iter = 1:100
    @time iterate(particle, isf, psi1, psi2, iter, omega, isJet)
end
