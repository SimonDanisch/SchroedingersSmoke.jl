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


isf = ISF(TorusDEC(vol_size..., vol_res...), hbar, dt)

# Set nozzle
const isJet = land(
    abs(isf.t.px - nozzle_cen[1]).<=nozzle_len/2,
    (isf.t.py - nozzle_cen[2]).^2+(isf.t.pz - nozzle_cen[3]).^2 .<= nozzle_rad.^2
)

# initialize psi
psi1 = ones(size(isf.t.px))
psi2 = psi1*0.01
psi1, psi2 = Normalize(psi1, psi2)

# constrain velocity
kvec = jet_velocity/isf.hbar;
omega = sum(jet_velocity.^2)/(2*isf.hbar);
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;

for iter = 1:10
    amp1 = abs(psi1)
    amp2 = abs(psi2)
    psi1 = 1.im * psi1
    psi2 = 1.im * psi2
    psi1[isJet] = amp1[isJet].*exp(1.im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1.im*phase[isJet])
    psi1, psi2 = PressureProject(isf, psi1, psi2)
end


## SET PARTICLES
particle = Particles(Float32[], Float32[], Float32[])


## MAIN ITERATION
itermax = ceil(tmax/dt);
function iterate(particle, isf, psi1, psi2, iter, omega, isJet)
    t = iter*dt
    # incompressible Schroedinger flow
    psi1, psi2 = SchroedingerFlow(isf, psi1,psi2)
    psi1, psi2 = Normalize(psi1,psi2)
    psi1, psi2 = PressureProject(isf, psi1,psi2)

    # constrain velocity
    phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz - omega*t;
    amp1 = abs(psi1);
    amp2 = abs(psi2);
    psi1[isJet] = amp1[isJet].*exp(1.im*phase[isJet])
    psi2[isJet] = amp2[isJet].*exp(1.im*phase[isJet])
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
end
#
using GLVisualize, GeometryTypes, GLWindow, GLAbstraction, Colors, GLFW
w=glscreen()
view(
    visualize(
        (Sphere(Point2f0(0), 0.005f0), (Float32[0], Float32[0], Float32[0])),
        color=RGBA{Float32}(0,0,0,0.3), billboard=true
    ),
    camera=:perspective
)

robj = renderlist(w)[1]

gpu_x, gpu_y, gpu_z = robj[:position_x], robj[:position_y], robj[:position_z]

for iter = 1:itermax
    isopen(w) || break
    @time iterate(particle, isf, psi1, psi2, iter, omega, isJet)
    update!(gpu_x, particle.x)
    update!(gpu_y, particle.y)
    update!(gpu_z, particle.z)
    GLWindow.render_frame(w)
    GLFW.PollEvents()
end
empty!(w)
yield()
GLFW.DestroyWindow(GLWindow.nativewindow(w))
