using SchroedingersSmoke

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
psi1 = (1.+0.0im)*psi1f
psi2 = (1.+0.0im)*psi2f
for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end


## SET PARTICLES
particle = Particles(Float32[], Float32[], Float32[])


using GLVisualize, GeometryTypes, GLWindow, GLAbstraction, Colors, GLFW
w=glscreen()
cubecamera(w)
# view(
#     visualize(
#         reinterpret(Vec3f0, isf.t.velocity),
#         ranges=map(x->0:x, vol_size),
#         color_norm = Vec2f0(0, 6),
#         color_map = RGBA{Float32}[RGBA{Float32}(0,1,0,0.6), RGBA{Float32}(1,0,0,1)]
#     ),
#     camera=:perspective
# )
empty!(w)
view(
    visualize(
        (Circle(Point2f0(0), 0.006f0), (Float32[0], Float32[0], Float32[0])),
        billboard=true
    ),
    camera=:perspective
)

#gpu_velocity = renderlist(w)[1][:rotation]
gpu_position_x = renderlist(w)[1][:position_x]
gpu_position_y = renderlist(w)[1][:position_y]
gpu_position_z = renderlist(w)[1][:position_z]


for iter = 1:2000
    isopen(w) || break
    t = iter*dt
    # incompressible Schroedinger flow
    psi1, psi2 = SchroedingerFlow(isf, psi1,psi2)
    psi1, psi2 = Normalize(psi1,psi2)
    psi1, psi2 = PressureProject(isf, psi1,psi2)

    # constrain velocity
    phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz - omega*t;
    amp1 = abs(psi1);
    amp2 = abs(psi2);
    psi1[isJet] = amp1[isJet].*exp(1.0im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1.0im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)

    # particle birth
    rt = rand(Float32, n_particles)*2*pi;
    newx = nozzle_cen[1]*ones(Float32, size(rt))
    newy = nozzle_cen[2] + 0.9*nozzle_rad*cos(rt)
    newz = nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)
    append!(particle.x, newx)
    append!(particle.y, newy)
    append!(particle.z, newz)
    # advect and show particles
    vx,vy,vz = VelocityOneForm(isf, psi1, psi2, isf.hbar);
    vx,vy,vz = StaggeredSharp(isf.t,vx,vy,vz);
    StaggeredAdvect(particle, isf.t,vx,vy,vz,isf.dt);
    Keep(particle,
        (particle.x .> 0f0) & (particle.x .< vol_size[1]) &
        (particle.y .> 0f0) & (particle.y .< vol_size[2]) &
        (particle.z .> 0f0) & (particle.z .< vol_size[3])
    )
    #update!(gpu_velocity, vec(reinterpret(Vec3f0, velocity)))
    update!(gpu_position_x, particle.x)
    update!(gpu_position_y, particle.y)
    update!(gpu_position_z, particle.z)
    render_frame(w)
    GLFW.PollEvents()

end
# empty!(w)
# yield()
# GLFW.DestroyWindow(GLWindow.nativewindow(w))
# create_video(frames, "test2", pwd(), 1)
