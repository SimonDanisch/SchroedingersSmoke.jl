using GLVisualize, GeometryTypes, GLWindow, GLAbstraction, Colors, GLFW
w=glscreen()
# cubecamera(w)
vol_size = (4,2,2);   # box size
vol_res = (64,32,32); # volume resolution
hbar = 0.1f0;            # Planck constant
dt = 1/48;             # time step
tmax = 50
# empty!(w)
view(
    visualize(
        zeros(Vec3f0, vol_res),
        ranges=map(x->0:x, vol_size),
        color_norm = Vec2f0(0, 6),
        color_map = RGBA{Float32}[RGBA{Float32}(0,1,0,0.6), RGBA{Float32}(1,0,0,1)]
    ),
    camera=:perspective
)
view(
    visualize(
        (Circle(Point2f0(0), 0.01f0), (Float32[0], Float32[0], Float32[0])),
        billboard=true
    ),
    camera=:perspective
)
gpu_velocity = renderlist(w)[2][:rotation]
gpu_position_x = renderlist(w)[3][:position_x]
gpu_position_y = renderlist(w)[3][:position_y]
gpu_position_z = renderlist(w)[3][:position_z]


Base.FFTW.set_num_threads(4)
blas_set_num_threads(4)


using SchroedingersSmoke, GeometryTypes

jet_velocity = [1,0,0]; # jet velocity
nozzle_cen = [vol_size...]/2; # nozzle center
nozzle_len = 0.5;                   # nozzle length
nozzle_rad = 0.6;                   # nozzle radius
n_particles = 50;   # number of particles

## INITIALIZATION

isf = ISF(TorusDEC(vol_size, vol_res), hbar, dt)

# initialize psi
psi1f = ones(size(isf.t.px))
psi2f = psi1f*0.1
psi1f, psi2f = Normalize(psi1f, psi2f)
psi1 = SchroedingersSmoke.AddCircle(isf, psi1f, nozzle_cen, [-1, 0,0], 0.6, isf.t.dx*5)
psi2 = (1f0+0f0*im)*psi2f
psi1, psi2 = PressureProject(isf, psi1, psi2)

# constrain velocity
kvec = jet_velocity/isf.hbar;
omega = sum(jet_velocity.^2)/(2*isf.hbar);
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;

# convert to complex
#psi1 = (1f0+0f0*im)*psi1f
psi2 = (1f0+0f0*im)*psi2f


itermax = ceil(tmax/dt);
function vis(psi1,psi2)
    r = (abs(psi1)^2 - abs(psi2)^2)*0.5+0.5
    g = -2*imag(conj(psi1)*psi2)*0.5+0.5
    b = 2*real(conj(psi1)*psi2)*0.5+0.5
    RGBA{Float32}(b, g, r, 1)
end
psi1slate = reshape(psi1[:, 16, 1:32], (64, 32))
psi2slate = reshape(psi2[:, 16, 1:32], (64, 32))
colormap = map(vis, psi1slate, psi2slate)
#view(visualize(colormap, primitive=Rectangle{Float32}(0,0, 1024, 1024)))
#colorgpu = w.renderlist[1][1][:image]
breakloop = false
particle = Particles(Float32[], Float32[], Float32[])
rt = rand(Float32, 1000)*2*pi;
newx = nozzle_cen[1]*ones(Float32, size(rt))
newy = nozzle_cen[2] + 0.9*nozzle_rad*cos(rt)
newz = nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)
append!(particle.x, newx)
append!(particle.y, newy)
append!(particle.z, newz)

for iter = 1:200
    println(iter)
    isopen(w) || breakloop || break
    #velocity = iterate(particle, isf, psi1, psi2, iter, omega)

    psi1, psi2 = SchroedingerFlow(isf, psi1, psi2)
    psi1, psi2 = Normalize(psi1,psi2)
    psi1, psi2 = PressureProject(isf, psi1,psi2)

    # particle birth

    # advect and show particles
    vx,vy,vz = VelocityOneForm(isf, psi1, psi2, isf.hbar);
    vx,vy,vz = StaggeredSharp(isf.t,vx,vy,vz);
    StaggeredAdvect(particle, isf.t,vx,vy,vz,isf.dt);

    Keep(particle,
        (particle.x .> 0f0) & (particle.x .< vol_size[1]) &
        (particle.y .> 0f0) & (particle.y .< vol_size[2]) &
        (particle.z .> 0f0) & (particle.z .< vol_size[3])
    )



    psi1slate = reshape(psi1[:, 16, 1:32], (64, 32))
    psi2slate = reshape(psi2[:, 16, 1:32], (64, 32))
    colormap = map(vis, psi1slate, psi2slate)

    update!(colorgpu, colormap)
    velocity = Vec3f0[Vec3f0(x...) for x in zip(vx,vy,vz)]
    update!(gpu_velocity, vec(reinterpret(Vec3f0, velocity)))
    update!(gpu_position_x, particle.x)
    update!(gpu_position_y, particle.y)
    update!(gpu_position_z, particle.z)
    
    render_frame(w)
    GLFW.PollEvents()
end

#@async renderloop(w)
velocity = iterate(particle, isf, psi1, psi2, 200, omega)
velocity = reshape(velocity, (64, 32,32))
x = reshape(velocity[:, 16, 1:32], (64, 32))
colormap = RGBA{Float32}[RGBA{Float32}(clamp((x[i,j]+1)./2, 0, 1)..., 1) for i=1:64, j=1:32]
