using GeometryTypes
using StaticArrays

function map_idx!{F, T}(f::F, A::AbstractArray, args::T)
    for x = 1:size(A, 1), y = 1:size(A, 2), z = 1:size(A, 3)
        idx = Vec{3, Int}(x, y, z)
        val = f(idx, A, args)
        @inbounds A[x, y, z] = val
    end
    A
end


@inline function inner_velocity_one_form(i, velocity, res_psi_hbar)
    res, psi1, psi2, hbar = res_psi_hbar
    i2 = mod.(i, res) + 1
    @inbounds begin
        psi12  = psi[i[1],  i[2],  i[3]]
        psix12 = psi[i2[1], i[2],  i[3]]
        psiy12 = psi[i[1],  i2[2] ,i[3]]
        psiz12 = psi[i[1],  i[2],  i2[3]]
    end
    psi1n = Vec(psix12[1], psiy12[1], psiz12[1])
    psi2n = Vec(psix12[2], psiy12[2], psiz12[2])
    angle.(
        conj(psi12[1]) .* psi1n .+
        conj(psi12[2]) .* psi2n
    ) * hbar
end
function velocity_one_form!(velocity, res, psi, hbar = 1.0f0)
    map_idx!(inner_velocity_one_form, velocity, (res, psi, hbar))
end



function gauge_transform!(psi, q)
    broadcast!(psi, psi, q) do psi, q
        eiq = exp(1.0im * (-q))
        (psi[1] * eiq, psi[2] * eiq)
    end
end

function div!(f, res, d, velocity)
    ds = inv.(d .^ 2)
    map_idx!(f) do xyz
        x, y, z = xyz
        ix, iy, iz = mod.(xyz - 2, res) + 1
        v1 = velocity[x,  y,  z]
        v2 = Vec(
            velocity[ix, y,  z][1],
            velocity[x,  iy, z][2],
            velocity[x,  y, iz][3]
        )
        sum((v1 .- v2) .* ds)
    end
    f
end


function poisson_solve(fac, f)
    fc = fft(f)
    fc .= (*).(fc, fac)
    ifft!(fc)
    fc
end

function pressure_project!(velocity, psi, fac, f, res, d)
    velocity_one_form!(velocity, res, psi)
    div!(f, res, d, velocity)
    q = poisson_solve(fac, f)
    gauge_transform!(psi, q)
end

function Normalize!(psi)
    broadcast!(psi, psi) do psi
        a, b = abs(psi[1]), abs(psi[2])
        norm = inv(sqrt(a*a + b*b))
        (psi[1] * norm, psi[2] * norm)
    end
    psi
end

function schroedinger_flow!(psi, mask)
    # extract single psi values into tmp arrays stored in obj
    psi1 = map(first, psi)
    psi2 = map(last, psi)
    psi1 = fftshift(fft(psi1)); psi2 = fftshift(fft(psi2));
    psi1 .= (*).(psi1, mask)
    psi2 .= (*).(psi2, mask)
    psi1 = ifft!(fftshift(psi1)); psi2 = ifft!(fftshift(psi2));
    psi .= tuple.(psi1, psi2) # convert back
    psi
end

type Particles
    xyz::Vector{Point3f0}
    length::Int
    active::Vector{Int}
end

function StaggeredAdvect(particle, velocity, dt, d, res, gridsize)
    map!(view(particle.xyz, particle.active)) do p

        k1 = staggered_velocity(velocity, p, d, gridsize, res)

        k2 = p + k1 .* dt ./ 2
        k2 = staggered_velocity(velocity, k2, d, gridsize, res)

        k3 = p + k2 .* dt ./ 2
        k3 = staggered_velocity(velocity, k3, d, gridsize, res)

        k4 = p + k3 .* dt
        k4 = staggered_velocity(velocity, k4, d, gridsize, res)

        p .+ dt/6 .* (k1 .+ 2*k2 .+ 2*k3 .+ k4)
    end
end

@inline function staggered_velocity(velocity, point, d, gs, res)
    p   = mod.(Vec(point), gs)
    i   = Vec{3, Int}(floor.(p ./ d) + 1)
    ip  = mod.(i, res) + 1

    v0  = velocity[i[1], i[2], i[3]]

    pxp = velocity[ip[1], i[2], i[3]]
    pyp = velocity[i[1], ip[2], i[3]]
    pzp = velocity[i[1], i[2], ip[3]]

    vn = Vec3f0(
        velocity[i[1], ip[2], ip[3]][1],
        velocity[ip[1], i[2], ip[3]][2],
        velocity[ip[1], ip[2], i[3]][3]
    )
    pp  = Vec3f0(pyp[1], pxp[2], pxp[3])
    pp2 = Vec3f0(pzp[1], pzp[2], pyp[3])

    w   = p - (i - 1) .* d
    w1  = Vec3f0(w[3], w[3], w[2])
    w2  = Vec3f0(w[2], w[1], w[1])
    return Point3f0(
        (1 - w1) .* ((1 - w2) .* v0 + w2 .* pp) +
        w1 .* ((1 - w2) .* pp2 + w2 .* vn)
    )
end



function Base.append!(p::Particles, xyz)
    @assert length(xyz) + p.length < length(p.xyz)
    range = (p.length+1):(p.length+length(xyz))
    p.xyz[range] = xyz
    p.length = (p.length+length(xyz))
    append!(p.active, range)
    nothing
end
Base.length(p::Particles) = p.length


@inline function isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(pos[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((pos[2] - nozzle_cen[2])^2 +
    (pos[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end

dims = (64, 32, 32)
grid_res = Vec(dims)
hbar = 0.1f0; dt = 1f0/48f0;
grid_size = Vec(4, 2, 2)
res = Vec(grid_res)
d = Vec3f0(grid_size ./ grid_res)

jet_velocity = Vec3f0(1, 0, 0)
nozzle_cen = Vec3f0(2-1.7, 1-0.034, 1+0.066)
nozzle_len = 0.5;
nozzle_rad = 0.5;
n_particles = 50
ranges = Vec(ntuple(3) do i
    linspace(0, grid_size[i], dims[i])
end)


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
velocity = zeros(Vec3f0, dims)
psi = [(one(Complex64), one(Complex64) * 0.1f0)
for i=1:dims[1], j=1:dims[2], k=1:dims[3]]

using BenchmarkTools



@benchmark inner_velocity_one_form($(Vec(1,2,3)), $(velocity), $(res), $(psi), $(1f0))
function VelocityOneForm(ix,iy, iz, resx, resy, resz, psi1, psi2, hbar=1.0)
    ixp = mod(ix, resx) + 1;
    iyp = mod(iy, resy) + 1;
    izp = mod(iz, resz) + 1;
    vx = angle(
        conj(psi1).*view(psi1, ixp,:,:) +
        conj(psi2).*view(psi2, ixp,:,:)
    );
    vy = angle(
        conj(psi1).*view(psi1,:,iyp,:) +
        conj(psi2).*view(psi2,:,iyp,:)
    )
    vz = angle(
        conj(psi1).*view(psi1,:,:,izp) +
        conj(psi2).*view(psi2,:,:,izp)
    )
    vx = vx*hbar;
    vy = vy*hbar;
    vz = vz*hbar;
    vx,vy,vz
end
psi1 = map(first, psi)
psi2 = map(last, psi)
@time VelocityOneForm(1:dims[1],1:dims[2],1:dims[3], res..., psi1, psi2, 1.0)
@code_llvm map_idx!(inner_velocity_one_form, velocity, res, psi, 1f0)

Normalize!(psi)

# constrain velocity
kvec = jet_velocity ./ hbar
omega = sum(jet_velocity.^2) / (2*hbar)

for iter = 1:10
    map_idx!(psi) do xyz
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
    pressure_project!(velocity, psi, fac, f_tmp, res, d)
end
using GLVisualize, GLAbstraction

w = glscreen(); @async GLWindow.waiting_renderloop(w)


particle = Particles(
    zeros(Point3f0, 10_000),
    0, Int[]
)

newp = Point3f0[begin
    rt = rand()*2*pi
    Point3f0(
        nozzle_cen[1] - 0.1,
        nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
        nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
    )
end for i=1:n_particles]

append!(particle, newp)




filter!(particle.active) do i
    p = particle.xyz[i]
    all(i-> p[i] in 0:grid_size[i], 1:3)
end

particle_vis = visualize(
    (Circle(Point2f0(0), 0.005f0), particle.xyz),
    indices = particle.active
)
_view(particle_vis, camera=:perspective)

velocity_vis = visualize(
    velocity,
    ranges = Tuple(ranges)
)
#_view(velocity_vis, camera=:perspective)
function in_grid(i, particle = particle)
    p = particle.xyz[i]
    for i=1:3
        p[1] >= 0 && p[1] <= grid_size[i] || return false
    end
    true
end


for i=1:1000
    t = i*dt
    # incompressable schroedinger flow
    schroedinger_flow!(psi, mask)
    Normalize!(psi)
    pressure_project!(velocity, psi, fac, f_tmp, res, d)

    # constrain velocity
    map_idx!(psi) do xyz
        x, y, z = xyz
        _psi = psi[x, y, z]
        pos = getindex.(ranges, xyz)
        if isjet(pos, nozzle_cen, nozzle_len, nozzle_rad)
            phase = sum(kvec .* pos) - omega*t
            return map(_psi) do p
                abs(p) * exp(1im * phase)
            end
        end
        _psi
    end
    pressure_project!(velocity, psi, fac, f_tmp, res, d)
    set_arg!(velocity_vis, :rotation, vec(velocity))

    newp = Point3f0[begin
        rt = rand()*2*pi
        Point3f0(
            nozzle_cen[1],
            nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
            nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
        )
    end for i=1:n_particles]

    append!(particle, newp)
    # advect and show particles
    velocity_one_form!(velocity, res, psi)
    velocity .= (.*).(velocity, Scalar(inv.(d)))
    filter!(in_grid, particle.active)
    particle.length = length(particle.active)
    StaggeredAdvect(particle, velocity, dt, d, grid_res, grid_size)
    set_arg!(particle_vis, :position, particle.xyz)
    set_arg!(particle_vis, :indices, map(x-> Cuint(x-1), particle.active))
    yield()
end




function test(velocity, res, psi, hbar)
    @inbounds for i=1:
    i2 = mod.(i, res) + 1
    psi12  = psi[i[1],  i[2],  i[3]]
    psix12 = psi[i2[1], i[2],  i[3]]
    psiy12 = psi[i[1],  i2[2] ,i[3]]
    psiz12 = psi[i[1],  i[2],  i2[3]]
    psi1n = Vec(psix12[1], psiy12[1], psiz12[1])
    psi2n = Vec(psix12[2], psiy12[2], psiz12[2])
    angle.(
        conj(psi12[1]) .* psi1n .+
        conj(psi12[2]) .* psi2n
    ) * hbar
end
end
