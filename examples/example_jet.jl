Base.FFTW.set_num_threads(4)
blas_set_num_threads(4)
using SchroedingersSmoke, Reactive, GLPlot
using GLVisualize, GeometryTypes, GLWindow, GLAbstraction, Colors, GLFW, ModernGL
GLPlot.init()
# example_jet
# An example of incompressible Schroedinger flow producing a jet.
#
## PARAMETERS


function GeometryTypes.isinside(point, lower, upper)
    @inbounds for (p,l,u) in zip(point, lower, upper)
        (p > l && p < u) || return false
    end
    true
end
# view(
#     visualize(
#         zeros(Vec3f0, size(isf.t.px)),
#         ranges=map(x->0:x, vol_size),
#         color_norm = Vec2f0(0, 6),
#         color_map = RGBA{Float32}[RGBA{Float32}(0,1,0,0.6), RGBA{Float32}(1,0,0,1)]
#     ),
#     camera=:perspective
# )
# empty!(w)

function filter2ind!(f, result, collections...)
    for (i, elem) in enumerate(zip(collections...))
        if f(elem...)
            push!(result, i)
        end
    end
    result
end

function main()
    const vol_size = (4,2,2);   # box size
    const vol_res = (64,32,32); # volume resolution
    const hbar = 0.1;           # Planck constant
    const dt = 1/48;            # time step
    const tmax = 50;            # max time

    const jet_velocity = [1,0,0]; # jet velocity

    const nozzle_cen = [2-1.7, 1-0.034, 1+0.066]; # nozzle center
    const nozzle_len = 0.5;                   # nozzle length
    const nozzle_rad = 0.5;                   # nozzle radius

    const n_particles = 10;   # number of particles


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
    psi1 = (1.+0.0im)*psi1f
    psi2 = (1.+0.0im)*psi2f
    for iter = 1:10
        amp1 = abs(psi1)
        amp2 = abs(psi2)
        psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
        psi1, psi2 = PressureProject(isf, psi1, psi2)
    end
    const max_particles = 20_000
    const max_history = 19
    update = Signal(true)
    active_particles = foldp(Array(Int, 0), update) do v0, _
        k = 1
        while !done(v0, k)
            i,k = next(v0, k)
            if !isinside(
                    (particle.x[i], particle.y[i], particle.z[i]),
                    (0,0,0), vol_size
                )
                splice!(v0, k-1)
                k -= 1
            end
        end
        v0
    end

    active_lines = foldp(zeros(GLuint, 0), active_particles) do line_indices, indices
        resize!(line_indices, length(indices) * max_history)
        for (i, ind) in enumerate(indices)
            i0 = i-1
            istart = i0*max_history+1
            for j=istart:(istart+max_history-1)
                line_indices[j] = j-1
            end
        end
        line_indices
    end
    ## SET PARTICLES
    bb = AABB{Float32}(Vec3f0(0), Vec3f0(vol_size))
    particle = Particles(
        zeros(Float32, max_particles),
        zeros(Float32, max_particles),
        zeros(Float32, max_particles),
        0, value(active_particles)
    )
    particle_vis = glplot(
        (Circle(Point3f0(0), 0.006f0), (particle.x, particle.y, particle.z)),
        boundingbox=bb,
        billboard=true, indices=active_particles
    ).children[]

    lines_ram = fill(Point3f0(0), max_history, max_particles)
    lines_color_ram = fill(RGBA{Float32}(0,0,0,0), max_history, max_particles)
    cmap = RGBA{Float32}[RGBA{Float32}(0,1,0,0.1), RGBA{Float32}(1,0,0,1)]
    # ma
    lines3d = glplot(
        vec(lines_ram), :lines,
        color = vec(lines_color_ram),
        dims=Vec{2,UInt32}(max_history, max_particles),
        boundingbox=bb,
        indices=active_lines, thickness=0.3f0,
        preferred_camera=:perspective
    )


    #gpu_velocity = renderlist(w)[1][:rotation]
    gpu_position_x = particle_vis[:position_x]
    gpu_position_y = particle_vis[:position_y]
    gpu_position_z = particle_vis[:position_z]

    @async for iter = 1:2000
        t = iter*dt
        # incompressible Schroedinger flow
        psi1, psi2 = SchroedingerFlow(isf, psi1,psi2)
        psi1, psi2 = Normalize(psi1,psi2)
        psi1, psi2 = PressureProject(isf, psi1, psi2)

        # constrain velocity
        phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz - omega*t
        amp1 = abs(psi1);
        amp2 = abs(psi2);
        psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
        psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
        psi1, psi2  = PressureProject(isf, psi1, psi2)

        # particle birth
        rt   = rand(Float32, n_particles)*2*pi
        newx = nozzle_cen[1] * ones(Float32, size(rt))
        newy = nozzle_cen[2] + 0.9*nozzle_rad*cos(rt)
        newz = nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)

        append!(particle, newx, newy, newz)
        push!(update, true)
        # advect and show particles
        vx,vy,vz = VelocityOneForm(isf, psi1, psi2, isf.hbar)
        vx,vy,vz = StaggeredSharp(isf.t, vx,vy,vz)
        StaggeredAdvect(particle, isf.t, vx,vy,vz, isf.dt)

        # update!(gpu_velocity, vec(map(Vec3f0, vx,vy,vz)))
        r = 1:length(particle)
        gpu_position_x[r] = particle.x[r]
        gpu_position_y[r] = particle.y[r]
        gpu_position_z[r] = particle.z[r]

        for i=particle.active
            p = Point3f0(particle.x[i], particle.y[i], particle.z[i])
            lines_ram[:, i] = circshift(lines_ram[:, i], 1)
            lines_ram[1, i] = p
            lines_ram[2, i] = p


            px = mod(p[1], isf.t.sizex)
            py = mod(p[2], isf.t.sizey)
            pz = mod(p[3], isf.t.sizez)

            ix = floor(Int, px/isf.t.dx) + 1
            iy = floor(Int, py/isf.t.dy) + 1
            iz = floor(Int, pz/isf.t.dz) + 1
            len = norm(Vec3f0(vx[i], vy[i], vz[i]))

            lines_color_ram[:, i] = circshift(lines_color_ram[:, i], 1)
            lines_color_ram[1, i] = color_lookup(cmap, len, 0,0.5)
        end

        set_arg!(lines3d, :vertex, vec(lines_ram))
        set_arg!(lines3d, :color, vec(lines_color_ram))
        yield()
    end
    #create_video(frames, "test2", pwd(), 1)
end

main()



using GLVisualize, GeometryTypes;w=glscreen();@async renderloop(w)
using Colors, GLAbstraction
empty!(w)
N = 10
r = linspace(0, 800, N)
positions = Point2f0[(x, 0) for x in r]
lines = fill(Point2f0(0), N, N)
lines[1, :] = positions

lines_color = fill(RGBA{Float32}(1,0,0,0), N, N)

linesobj = visualize(
    lines, :lines, dims=size(lines),
    color = vec(lines_color), # needs to be 1D right now, to lessen the amount of automatic conversions
    boundingbox=nothing,
    thickness=3.0
)
_view(linesobj, camera=:orthographic_pixel)


for t in 0.1:0.1:0.9
    c = RGBA{Float32}(0, 0, t, 1)
    for i=1:length(positions)
        lines[:, i] = circshift(view(lines, :, i), 1)
        lines[1, i] = Point2f0(r[i], t*300)
        lines_color[:, i] = circshift(view(lines_color, :, i), 1)
        lines_color[1, i] = c
    end
    set_arg!(linesobj, :vertex, vec(lines))
    set_arg!(linesobj, :color, vec(lines_color))
    sleep(0.1)
    yield()
end
map(Tuple, lines[end, :])
