include("../src/MatlabPort.jl")
using GeometryTypes, MatlabPort
import MatlabPort: ISF, Normalize, PressureProject, land

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

isf = ISF(vol_size, vol_res, hbar, dt)

isJet = land(
    abs(isf.px - nozzle_cen[1]).<=nozzle_len/2,
    (isf.py - nozzle_cen[2]).^2+(isf.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
)

# initialize psi
psi1f = ones(size(isf.px))
psi2f = psi1f * 0.01
psi1f, psi2f = Normalize(psi1f, psi2f)

# constrain velocity
kvec = jet_velocity / isf.hbar
omega = sum(jet_velocity .^ 2) / (2*isf.hbar);
phase = kvec[1].*isf.px + kvec[2].*isf.py + kvec[3].*isf.pz;
psi1 = (1.+0.0im)*psi1f
psi2 = (1.+0.0im)*psi2f
for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet] .* exp(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet] .* exp(1.0im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end
