include("../src/MatlabPort.jl")
import MatlabPort

include("../src/BroadcastPort.jl")
import BroadcastPort
macro _t(name, x)
    quote
        tic()
        alloc = @allocated $(esc(x))
        t = toq()
        println($name, ": alloc: ", alloc, " t: ", t)
    end
end
using GeometryTypes, StaticArrays, Base.Test

vol_size = (4,2,2)   # box size
dims = (64, 32, 32) # volume resolution
hbar = 0.1           # Planck constant
dt = 1/48            # time step
isf_1 = MatlabPort.ISF(vol_size, dims, hbar, dt)
isf_2 = BroadcastPort.ISF{Int, Float32}(vol_size, dims, hbar, dt)
@testset "Schroedinger Mask" begin
    @test all(isapprox.(isf_1.mask, isf_2.mask))
end

psi1f = ones(dims)
psi2f = psi1f * 0.01
psi1f, psi2f = MatlabPort.Normalize(psi1f, psi2f)
psi1, psi2 = (1.+0.0im)*psi1f, (1.+0.0im)*psi2f

psi = [(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]]
BroadcastPort.Normalize!(psi)

@testset "initial condition" begin
    @test all(isapprox.(map(first, psi1), psi1))
    @test all(isapprox.(map(last, psi2), psi2))
end

jet_velocity = Vec3f0(1, 0, 0)
nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
nozzle_len = 0.5f0
nozzle_rad = 0.5f0

isJet = MatlabPort.land(
    abs(isf_1.px - nozzle_cen[1]).<=nozzle_len/2,
    (isf_1.py - nozzle_cen[2]).^2+(isf_1.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
)

kvec = jet_velocity / isf_1.hbar
omega = sum(jet_velocity .^ 2) / (2*isf_1.hbar);
phase = kvec[1].*isf_1.px + kvec[2].*isf_1.py + kvec[3].*isf_1.pz;
for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet] .* exp(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet] .* exp(1.0im*phase[isJet])
    psi1, psi2 = MatlabPort.PressureProject(isf_1, psi1, psi2)
end

function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end
kvec = jet_velocity ./ hbar
omega = sum(jet_velocity.^2f0) / (2f0*hbar)

# constrain velocity
for iter = 1:10
    BroadcastPort.map_idx!(psi, ()) do xyz, psi, _
        x, y, z = xyz
        psie = psi[x, y, z]
        pos = getindex.(isf_2.ranges, xyz)
        if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
            phase = sum(kvec .* pos)
            return map(psie) do p
                abs(p) * exp(1im * phase)
            end
        end
        psie
    end
    BroadcastPort.pressure_project!(isf_2, psi)
end
@testset "constrain" begin
    @test all(isapprox.(map(first, psi1), psi1))
    @test all(isapprox.(map(last, psi2), psi2))
end

for i=1:10
    psi1, psi2 = MatlabPort.SchroedingerFlow(isf_1, psi1, psi2)
    psi1, psi2 = MatlabPort.Normalize(psi1, psi2)
    psi1, psi2 = MatlabPort.PressureProject(isf_1, psi1, psi2)
end

for i=1:10
    BroadcastPort.schroedinger_flow!(isf_2, psi)
    BroadcastPort.Normalize!(psi)
    BroadcastPort.pressure_project!(isf_2, psi)
end
@testset "pumping the flow" begin
    @test all(isapprox.(map(first, psi1), psi1))
    @test all(isapprox.(map(last, psi2), psi2))
end

N = 200
particle_1 = MatlabPort.Particles(
    rand(N)*vol_size[1],
    rand(N)*vol_size[2],
    rand(N)*vol_size[3],
    N, collect(1:N)
)
particles = map(Point3f0, zip(particle_1.x, particle_1.y, particle_1.z))

# this should not change anything but give us vx,vy,vz
vx,vy,vz = MatlabPort.VelocityOneForm(isf_1, psi1, psi2, isf_1.hbar)

for i=1:1
    @testset "simulation loop i = $i" begin
        t = i * dt
        psi1, psi2 = MatlabPort.SchroedingerFlow(isf_1, psi1, psi2)
        psi1, psi2 = MatlabPort.Normalize(psi1, psi2)
        psi1, psi2 = MatlabPort.PressureProject(isf_1, psi1, psi2)

        @_t "flow" BroadcastPort.schroedinger_flow!(isf_2, psi)
        @_t "normalize" BroadcastPort.Normalize!(psi)
        @_t "pressure prj" BroadcastPort.pressure_project!(isf_2, psi)

        @test all(isapprox.(map(first, psi1), psi1))
        @test all(isapprox.(map(last, psi2), psi2))

        # constrain velocity

        phase = kvec[1].*isf_1.px + kvec[2].*isf_1.py + kvec[3].*isf_1.pz - omega*t
        amp1 = abs(psi1)
        amp2 = abs(psi2)
        psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
        psi1, psi2  = MatlabPort.PressureProject(isf_1, psi1, psi2)

        @_t "restrict" BroadcastPort.map_idx!(psi, ()) do xyz, psi, _
            x, y, z = xyz
            _psi = psi[x, y, z]
            pos = getindex.(isf_2.ranges, xyz)
            if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
                phase = sum(kvec .* pos) - omega*t
                return map(_psi) do p
                    abs(p) * exp(1im * phase)
                end
            end
            _psi
        end
        @_t "pressure prj" BroadcastPort.pressure_project!(isf_2, psi)

        @test all(isapprox.(map(first, psi1), psi1))
        @test all(isapprox.(map(last, psi2), psi2))

        vx,vy,vz = MatlabPort.VelocityOneForm(isf_1, psi1, psi2, isf_1.hbar)
        @_t "veloneform" BroadcastPort.velocity_one_form!(isf_2, psi, isf_2.hbar)
        @show mean(map(Vec3f0, zip(vx, vy, vz)) - isf_2.velocity)
        @testset "velocity after one form" begin
            @test all(isapprox.(map(Vec3f0, zip(vx, vy, vz)), isf_2.velocity))
        end
        vx,vy,vz = MatlabPort.StaggeredSharp(isf_1, vx,vy,vz)
        @_t "stagger" map!(v-> v ./ isf_2.d, isf_2.velocity) # inplace StaggeredSharp
        @testset "velocity after staggered sharp" begin
            @test all(isapprox.(map(Vec3f0, zip(vx, vy, vz)), isf_2.velocity))
        end

        MatlabPort.StaggeredAdvect(particle_1, isf_1, vx,vy,vz, isf_1.dt)
        @_t "advect" BroadcastPort.staggered_advect!(particles, isf_2)

        particlep = map(Point3f0, zip(particle_1.x, particle_1.y, particle_1.z))
        @test all(isapprox.(particles, particlep))
    end
end
