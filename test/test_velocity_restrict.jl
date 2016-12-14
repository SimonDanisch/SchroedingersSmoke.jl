include("../src/MatlabPort.jl")
using GeometryTypes, MatlabPort
import MatlabPort: ISF, Normalize, PressureProject, land
using GeometryTypes
using StaticArrays

function test_matport()
    vol_size = (4,2,2)   # box size
    vol_res = (64,32,32) # volume resolution
    hbar = 0.1           # Planck constant
    dt = 1/48            # time step
    jet_velocity = [1,0,0]; # jet velocity
    nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; # nozzle center
    nozzle_len = 0.5;                   # nozzle length
    nozzle_rad = 0.5;                   # nozzle radius

    isf = MatlabPort.ISF(vol_size, vol_res, hbar, dt)

    isJet = MatlabPort.land(
        abs(isf.px - nozzle_cen[1]).<=nozzle_len/2,
        (isf.py - nozzle_cen[2]).^2+(isf.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
    )

    # initialize psi
    psi1f = ones(size(isf.px))
    psi2f = psi1f * 0.01
    psi1f, psi2f = MatlabPort.Normalize(psi1f, psi2f)

    # constrain velocity
    kvec = jet_velocity / isf.hbar
    phase = kvec[1].*isf.px + kvec[2].*isf.py + kvec[3].*isf.pz;
    psi1 = (1.+0.0im)*psi1f
    psi2 = (1.+0.0im)*psi2f
    for iter = 1:10
        amp1 = abs(psi1)
        amp2 = abs(psi2)
        psi1[isJet] = amp1[isJet] .* exp(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet] .* exp(1.0im*phase[isJet])
        psi1, psi2 = MatlabPort.PressureProject(isf, psi1, psi2)
    end
    psi1, psi2
end


function isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(pos[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((pos[2] - nozzle_cen[2])^2 +
    (pos[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end
# testing the generic broadcast implementation
include("../src/BroadcastPort.jl")
import BroadcastPort
function test_broadcast()
    dims = (64, 32, 32)
    grid_res = Vec(dims)
    hbar = 0.1f0; dt = 1f0/70f0;
    grid_size = Vec(4, 2, 2)
    res = Vec(grid_res)
    d = Vec3f0(grid_size ./ grid_res)
    ranges = Vec(ntuple(3) do i
        ((1:grid_res[i]) - 1) * d[i]
    end)
    idx2 = ntuple(3) do i
        mod(1:res[i], res[i]) + 1
    end
    psi = [(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]
    BroadcastPort.Normalize!(psi)
    jet_velocity = Vec3f0(1, 0, 0)
    nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
    nozzle_len = 0.5;
    nozzle_rad = 0.5;
    kvec = jet_velocity ./ hbar
    omega = sum(jet_velocity.^2) / (2*hbar)
    velocity = zeros(Vec3f0, dims)
    fac = zeros(Float32, dims)
    mask = zeros(Complex64, dims)
    f = -4 * pi^2 * hbar
    for x = 1:dims[1], y = 1:dims[2], z = 1:dims[3]
        xyz = Vec(x, y, z)
        # fac
        s = sin(pi * (xyz - 1) ./ res) ./ d
        denom = sum(s .^ 2)
        fac[x,y,z] = -0.25f0 ./ denom

        # schroedingers mask
        k = (xyz - 1 - res ./ 2) ./ grid_size
        lambda = f * sum(k .^ 2)
        mask[x,y,z] = exp(1.0im * lambda * dt / 2)
    end
    fac[1,1,1] = 0

    f_tmp = zeros(Float32, dims)
    for iter = 1:10
        BroadcastPort.map_idx!(psi, ()) do xyz, psi, _
            x, y, z = xyz
            _psi = psi[x, y, z]
            pos = getindex.(ranges, xyz)
            if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
                phase = sum(kvec .* pos)
                return map(_psi) do p
                    abs(p) * exp(1im * phase)
                end
            end
            _psi
        end
        BroadcastPort.pressure_project!(velocity, psi, idx2, fac, f_tmp, res, d)
    end
    map(first, psi), map(last, psi)
end
apsi1, apsi2 = test_matport()

bpsi1, bpsi2 = test_broadcast()


vol_size = (4,2,2)   # box size
vol_res = (64,32,32) # volume resolution
hbar = 0.1           # Planck constant
dt = 1/48            # time step
jet_velocity = [1,0,0]; # jet velocity
nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; # nozzle center
nozzle_len = 0.5;                   # nozzle length
nozzle_rad = 0.5;                   # nozzle radius
dims = vol_res
isf = MatlabPort.ISF(vol_size, vol_res, hbar, dt)
d = Vec3f0(vol_size) ./ Vec3f0(vol_res)

for i=1:dims[1], j=1:dims[2], k=1:dims[3]
    isf.px[i, j, k] == ranges[1][i] &&
    isf.py[i, j, k] == ranges[2][j] &&
    isf.pz[i, j, k] == ranges[3][k] || error("lol")
end
