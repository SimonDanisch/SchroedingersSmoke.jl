Base.FFTW.set_num_threads(8)

using SchroedingersSmoke
import SchroedingersSmoke.ParallelPort
using StaticArrays, Colors
x = code_typed(map, Tuple{typeof(inv), StaticArrays.SVector{3,Float32}}, optimize = false)

import ParallelPort: ISF, normalize_psi, pressure_project!
import ParallelPort: velocity_one_form!, schroedinger_flow!
import ParallelPort: Particles, staggered_advect!, map_idx!
import ParallelPort: Vec, Vec3f0, Point, Point3f0
x = code_typed(map, Tuple{typeof(inv), StaticArrays.SVector{3, Float32}}, optimize = false)

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
        return map(psi) do p
            abs(p) * exp(1im * phase)
        end
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
psi = [(one(Complex64), one(Complex64) * 0.01f0) for i=1:dims[1], j=1:dims[2], k=1:dims[3]];
normalize_psi.(psi);

kvec = jet_velocity ./ hbar;
omega = sum(jet_velocity.^2f0) / (2f0*hbar);

# constrain velocity
for iter = 1:10
    restrict_velocity!(isf, psi, kvec, nozzle_cen, nozzle_len, nozzle_rad)
    pressure_project!(isf, psi)
end
using GPUArrays
using GPUArrays.JLBackend
using GPUArrays.JLBackend.JLArray

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
#append!(particle, newp);

function simloop(
        N, isf, psi, kvec, omega, n_particles,
        nozzle_rad, nozzle_cen, particle
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
        # newp = map(1:n_particles) do _
        #     rt = rand()*2*pi
        #     Point3f0(
        #         nozzle_cen[1] - 0.1,
        #         nozzle_cen[2] + 0.9 * nozzle_rad * cos(rt),
        #         nozzle_cen[3] + 0.9 * nozzle_rad * sin(rt)
        #     )
        # end
        #append!(particle, newp)
        #filter!(in_grid, particle.active)

        # advect and show particles
        velocity_one_form!(isf, psi, isf.hbar)
        # inplace StaggeredSharp
        dinv = inv.(isf.d)
        broadcast!(x-> x .* dinv, isf.velocity, isf.velocity)

        staggered_advect!(view(particle.xyz, particle.active), isf)
    end
end
using Sugar
import Sugar: @lazymethod, getsource!, dependencies!, getast!, isfunction, istype

import Sugar: LazyMethod

immutable ResolveError <: Exception
    stack::Vector
end
function Base.showerror(io::IO, e::ResolveError)
    println(io, "$(e.stack[1])\nIn:")
    for elem in e.stack[2:end]
        println(io, "  ", elem)
    end
end

function resolve_dependencies!(dep::LazyMethod, visited = LazyMethod(Void), stack = [])
    Sugar.isintrinsic(dep) && return visited.dependencies
    if isa(dep.signature, Module)
        delete!(visited.dependencies, dep)
        return visited.dependencies
    end
    if dep in visited.dependencies
        # when already in deps we need to move it up!
        delete!(visited.dependencies, dep)
        push!(visited.dependencies, dep)
    else
        push!(visited, dep)
        try
            push!(stack, dep.signature)
            resolve_dependencies!(dependencies!(dep), visited)

        catch e
            for elem in stack
                println("  ", elem)
            end
        finally
            pop!(stack)
        end
    end
    visited.dependencies
end
function resolve_dependencies!(deps, visited, stack = [])
    for dep in copy(deps)
        if !Sugar.isintrinsic(dep)
            push!(stack, dep.signature)
            resolve_dependencies!(dep, visited)
            pop!(stack)
        end
    end
    visited.dependencies
end
x = @lazymethod simloop(
   100, isf, psi, kvec, omega, n_particles,
   nozzle_rad, nozzle_cen, particle
)
dependencies!(x)

deps = resolve_dependencies!(x);
length(deps)
schroefun_all = filter(deps) do fun
    try
        return (Sugar.getmethod(fun).module == SchroedingersSmoke)
    end
    false
end
typs = filter(Sugar.istype, deps)
abst_t = filter(typs) do t
    !isleaftype(t.signature)
end
println(length(abst_t))
for elem in schroefun_all
    try
        println(elem.signature)
        typs = filter(Sugar.istype, dependencies!(elem))
        abst_t = collect(filter(typs) do t
            !isleaftype(t.signature)
        end)
        println(length(abst_t))
    end
end
x = code_typed(map, Tuple{typeof(inv), StaticArrays.SVector{3,Float32}}, optimize = false)

lol = LazyMethod((
    Base._similar_for,
    Tuple{
        Array{Tuple{Complex{Float32},Complex{Float32}},3},
        Type{Complex{Float32}},
        Base.Generator{Array{Tuple{Complex{Float32},Complex{Float32}},3}, typeof(first)},
        Base.HasShape
    }
))
ast = getast!(lol)
heh = dependencies!(lol)[4]



x = dependencies!(lol)[4]
lal = LazyMethod((Base.similar, Tuple{
        Array{Tuple{Complex{Float32},Complex{Float32}}, 3},
        Int64,
        Tuple{Int64,Int64,Int64}
    }
))
getast!(lal)

lal = LazyMethod((similar,
    Tuple{
        Array{Tuple{Complex{Float32},Complex{Float32}},3},
        DataType,
        Tuple{Int64,Int64,Int64}
    }
))
ast = Sugar.get_ast(code_typed, similar,
    Tuple{
        Array{Tuple{Complex{Float32},Complex{Float32}},3},
        DataType,
        Tuple{Int64,Int64,Int64}
    }
)

# args = (100, isf, psi, kvec, omega, n_particles,
# nozzle_rad, nozzle_cen, particle)
# targs = map(typeof, args)
# ast = Sugar.get_ast(code_typed, simloop, targs)
# using MacroTools
# MacroTools.walk
# str = sprint() do io
#     show(io, Expr(:block, ast...))
# end
# test = nothing
# for elem in ast
#     if isa(elem, Expr) && elem.head == :(=) && SlotNumber(17) == elem.args[1]
#         println(typeof(elem.args[1]))
#         println("lel ", elem)
#         test = elem
#     end
# end
# test.args[2].args[1]
