type ISF{T}
    t::TorusDEC
    hbar::T             # reduced Planck constant
    dt::T               # time step
    SchroedingerMask::cl.CLArray{Complex{T},3} # Fourier coefficient for solving Schroedinger eq
    psi::NTuple{2, cl.CLArray{Complex{T}, 3}}
    velocity::cl.CLArray{Vec3f0, 3}
    function ISF(t, hbar, dt)
        isf = new{T}(t, hbar, dt)
        BuildSchroedinger(isf)
        isf
    end
end
ISF{T}(t, hbar::T, dt::T) = ISF{T}(t, hbar, dt)

"""
builds coefficients in Fourier space.
"""
function BuildSchroedinger(obj::ISF)
    nx=obj.t.resx; ny=obj.t.resy; nz=obj.t.resz
    fac = -4*pi^2*obj.hbar
    kx = (obj.t.iix-1-nx/2)/(obj.t.sizex)
    ky = (obj.t.iiy-1-ny/2)/(obj.t.sizey)
    kz = (obj.t.iiz-1-nz/2)/(obj.t.sizez)
    lambda = fac*(kx.^2+ky.^2+kz.^2)
    obj.SchroedingerMask = CLArray(ifftshift(exp(1f0*im*lambda*obj.dt/2.)))
end

"""
solves Schroedinger equation for dt time.
"""
function SchroedingerFlow(obj, psi1, psi2)
    psi1 = fft!(psi1)
    psi2 = fft!(psi2)
    mul!(psi1, obj.SchroedingerMask)
    mul!(psi2, obj.SchroedingerMask)
    ifft!(psi1)
    ifft!(psi2)
    nothing
end

"""
Pressure projection of 2-component wave function.
"""
function PressureProject(obj, psi1, psi2)
    velocity = VelocityOneForm(obj, psi1,psi2)
    div = Div(obj.t, velocity)
    q = PoissonSolve(obj.t, div)
    GaugeTransform(psi1, psi2, -q)
end


@cl_kernel VelocityOneForm program_ISF psi1 psi2 velocity hbar

function VelocityOneForm(obj, psi1, psi2, hbar=1.0f0)
    velocity = obj.t.velocity
    psi1, psi2 = obj.psi
    VelocityOneForm(psi1, psi2, velocity, hbar)
end


@cl_kernel GaugeTransform program_ISF psi1 psi2 q


@cl_kernel Normalize program_ISF psi1 psi2
