using SchroedingersSmoke, GeometryTypes

Base.FFTW.set_num_threads(4)
Base.BLAS.set_num_threads(4)

vol_size = (4,2,2);   # box size
vol_res = (64,32,32); # volume resolution
hbar = 0.1f0;         # Planck constant
dt = 1f0/48f0;        # time step
jet_velocity = [1,0,0]; # jet velocity
nozzle_cen = [vol_size...]/2; # nozzle center
nozzle_len = 0.5;                   # nozzle length
nozzle_rad = 0.6;                   # nozzle radius
n_particles = 10;   # number of particles
const max_particles = 20_000

## INITIALIZATION

isf = ISF(TorusDEC(vol_size, vol_res), hbar, dt)

# initialize psi
psi1f = ones(size(isf.t.px))
psi2f = psi1f*0.1
psi1f, psi2f = Normalize(psi1f, psi2f)
psi2 = (1f0+0f0*im)*psi2f
psi1, psi2 = PressureProject(isf, psi1f, psi2f)

# constrain velocity
kvec = jet_velocity/isf.hbar;
omega = sum(jet_velocity.^2)/(2*isf.hbar);
phase = kvec[1].*isf.t.px + kvec[2].*isf.t.py + kvec[3].*isf.t.pz;

#

# convert to complex
#psi1 = (1f0+0f0*im)*psi1f
psi2 = (1f0+0f0*im)*psi2f

itermax = ceil(200/dt);
breakloop = false

particle = Particles(
    zeros(Point3f0, max_particles),
    0, [1,2,3,4,5]
)

newp = Point3f0[begin
    rt = rand()*2*pi
    Point3f0(
        nozzle_cen[1],
        nozzle_cen[2] + 0.9*nozzle_rad*cos(rt),
        nozzle_cen[3] + 0.9*nozzle_rad*sin(rt)
    )
end for i=1:n_particles]

append!(particle.xyz, newp)

for iter = 1:10
    psi1, psi2 = SchroedingerFlow(isf, psi1, psi2)
    psi1, psi2 = Normalize(psi1,psi2)
    psi1, psi2 = PressureProject(isf, psi1,psi2)

    # advect and show particles
    velocity = VelocityOneForm(isf, psi1, psi2, isf.hbar);
    velocity = StaggeredSharp(isf.t,velocity);
    StaggeredAdvect(particle, isf.t, velocity,isf.dt);

end
