# ISF a handle class for simulating incompressible Schroedinger flow.
#
# SYNTAX
#
#   isf = ISF(sizex,sizey,sizez,resx,resy,resz)
#   isf = ISF(sizex,sizey,sizez,res)
#
# DESCRIPTION
#
#   ISF is a subclass of TorusDEC handle class. See TorusDEC for calling
#   constructor.
#
#   To setup,
#
#     isf = ISF(5,5,5,32,32,32);  # call constructor
#     isf.hbar = 0.1;             # specify Planck constant
#     isf.dt   = 1/24;            # specify time step
#     isf.BuildSchroedinger;      # this command builds coeff for solving
#                                 # Schroedinger equation in Fourier domain
#
#   Useful functions:
#
#     [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2)
#         solves Schroedinger equation for (psi1,psi2) for isf.dt time.
#
#     [psi1,psi2] = isf.Normalize(psi1,psi2)
#         normalizes (psi1,psi2).
#
#     [psi1,psi2] = isf.PressureProject(psi1,psi2)
#         projects (psi1,psi2) to satisfy divergence free constraint.
#
#     [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2,isf.hbar)
#         extracts velocity 1-from from (psi1,psi2)
#
#     [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2)
#         extracts velocity 1-form assuming hbar=1.
#
#
# See also TorusDEC
type ISF{T}
    t::TorusDEC
    hbar::T             # reduced Planck constant
    dt::T               # time step
    SchroedingerMask::Array{Complex{T},3} # Fourier coefficient for solving Schroedinger eq
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
    obj.SchroedingerMask = exp(1.0im*lambda*obj.dt/2.)
end

"""
solves Schroedinger equation for dt time.
"""
function SchroedingerFlow(obj, psi1, psi2)
    psi1 = fftshift(fft(psi1)); psi2 = fftshift(fft(psi2));
    @inbounds for i=1:length(psi1)
        psi1[i] = psi1[i]*obj.SchroedingerMask[i];
        psi2[i] = psi2[i]*obj.SchroedingerMask[i];
    end
    psi1 = ifft!(fftshift(psi1)); psi2 = ifft!(fftshift(psi2));
    psi1, psi2
end

"""
Pressure projection of 2-component wave function.
"""
function PressureProject(obj, psi1, psi2)
    velocity = VelocityOneForm(obj, psi1,psi2)
    div = Div(obj.t, velocity)
    q = PoissonSolve(obj.t, div);
    GaugeTransform(psi1, psi2, -q)
end

"""
extracts velocity 1-form from (psi1,psi2).
If hbar argument is empty, hbar=1 is assumed.
"""
function VelocityOneForm(obj, psi1, psi2, hbar=1.0)
    @inbounds for z=obj.t.iz, y=obj.t.iy
        @simd for x=obj.t.ix
            ixp = mod(x, obj.t.resx) + 1
            iyp = mod(y, obj.t.resy) + 1
            izp = mod(z, obj.t.resz) + 1
            psi1c = conj(psi1[x,y,z]); psi2c = conj(psi1[x,y,z]);
            vx = angle(psi1c*psi1[ixp,y,z] + psi2c*psi2[ixp,y,z])
            vy = angle(psi1c.*psi1[x,iyp,z] + psi2c.*psi2[x,iyp,z])
            vz = angle(psi1c.*psi1[x,y,izp] + psi2c.*psi2[x,y,izp])
            obj.t.velocity[x,y,z] = Point3f0(vx*hbar, vy*hbar, vz*hbar)
        end
    end
    obj.t.velocity
end


"""
adds a vortex ring to a 1-component wave function psi.
Inputs center, normal, r specify the circle.
Input d specify the thickness around the disk to create a boost
in phase. Usually d = 5*dx where dx is grid edge length.
"""
function AddCircle(obj, psi, center, normal, r, d)
    rx = obj.t.px - center[1]
    ry = obj.t.py - center[2]
    rz = obj.t.pz - center[3]
    normal = normal/norm(normal,2);
    alpha = zeros(size(rx));
    z = rx*normal(1) + ry*normal(2) + rz*normal(3);
    inCylinder = rx.^2+ry.^2+rz.^2 - z.^2 < r^2;
    inLayerP = z> 0 & z<= d/2 & inCylinder
    inLayerM = z<=0 & z>=-d/2 & inCylinder
    alpha[inLayerP] = -pi*(2*z(inLayerP)/d - 1)
    alpha[inLayerM] = -pi*(2*z(inLayerM)/d + 1)
    psi.*exp(1.im*alpha);
end

"""
multiplies exp(i*q) to (psi1,psi2)
"""
function GaugeTransform(psi1, psi2, q)
    @inbounds for i=1:length(psi1)
        psi1[i], psi2[i] = _gauge_transform(psi1[i], psi2[i], q[i])
    end
    psi1, psi2
end
function _gauge_transform(psi1, psi2, q)
    eiq = exp(1.0im*q);
    psi1 = psi1*eiq
    psi2 = psi2*eiq
    psi1,psi2
end

"""
extracts Clebsch variable s=(sx,sy,sz) from (psi1,psi2)
"""
function Hopf(psi1,psi2)
    nd = size(psi1)
    sx = Array(Float64, nd)
    sy = Array(Float64, nd)
    sz = Array(Float64, nd)
    Hopf(psi1,psi2,sx,sy,sz)
end
function Hopf(psi1,psi2,sx,sy,sz)
    @inbounds for i=1:length(sx)
        sx[i], sy[j], sz[j] = _hopf(psi1[i], psi2[i])
    end
    sx,sy,sz
end
function _hopf(psi1,psi2)
    a = real(psi1)
    b = imag(psi1)
    c = real(psi2)
    d = imag(psi2)
    sx = 2*(a*c + b*d)
    sy = 2*(a*d - b*c)
    sz = a^2 + b^2 - c^2 - d^2
    sx,sy,sz
end

Base.@pure function _normalize(psi1, psi2)
    norm = sqrt(abs(psi1)^2 + abs(psi2)^2)
    psi1/norm, psi2/norm
end
"""
normalizes (psi1,psi2)
"""
function Normalize(psi1,psi2)
    @inbounds for (i,j) in zip(eachindex(psi1), eachindex(psi2))
        psi1[i], psi2[j] = _normalize(psi1[i], psi2[j])
    end
    psi1, psi2
end
