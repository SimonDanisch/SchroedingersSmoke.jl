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
type ISF
    t::TorusDEC
    hbar::Float64             # reduced Planck constant
    dt::Float64               # time step
    SchroedingerMask::Array{Complex{Float64},3} # Fourier coefficient for solving Schroedinger eq
    function ISF(t, hbar, dt)
        isf = new(t, hbar, dt)
        BuildSchroedinger(isf)
        isf
    end
end

"""
builds coefficients in Fourier space.
"""
function BuildSchroedinger(obj::ISF)
    nx=obj.t.resx; ny=obj.t.resy; nz=obj.t.resz;
    fac = -4*pi^2*obj.hbar;
    kx = (obj.t.iix-1-nx/2)/(obj.t.sizex);
    ky = (obj.t.iiy-1-ny/2)/(obj.t.sizey);
    kz = (obj.t.iiz-1-nz/2)/(obj.t.sizez);
    lambda = fac*(kx.^2+ky.^2+kz.^2);
    obj.SchroedingerMask = exp(1.0im*lambda*obj.dt/2.);
end

"""
solves Schroedinger equation for dt time.
"""
function SchroedingerFlow(obj, psi1, psi2)
    psi1 = fftshift(fft(psi1)); psi2 = fftshift(fft(psi2));
    psi1 = psi1.*obj.SchroedingerMask;
    psi2 = psi2.*obj.SchroedingerMask;
    psi1 = ifft(fftshift(psi1)); psi2 = ifft(fftshift(psi2));
    psi1, psi2
end

"""
Pressure projection of 2-component wave function.
"""
function PressureProject(obj, psi1, psi2)
    vx,vy,vz = VelocityOneForm(obj, psi1,psi2)
    div = Div(obj.t, vx,vy,vz)
    q = PoissonSolve(obj.t, div);
    GaugeTransform(psi1, psi2, -q)
end

"""
extracts velocity 1-form from (psi1,psi2).
If hbar argument is empty, hbar=1 is assumed.
"""
function VelocityOneForm(obj, psi1, psi2, hbar=1.0)
    ixp = mod(obj.t.ix, obj.t.resx) + 1;
    iyp = mod(obj.t.iy, obj.t.resy) + 1;
    izp = mod(obj.t.iz, obj.t.resz) + 1;
    vx = angle(
        conj(psi1).*sub(psi1, ixp,:,:) +
        conj(psi2).*sub(psi2, ixp,:,:)
    );
    vy = angle(
        conj(psi1).*sub(psi1,:,iyp,:) +
        conj(psi2).*sub(psi2,:,iyp,:)
    )
    vz = angle(
        conj(psi1).*sub(psi1,:,:,izp) +
        conj(psi2).*sub(psi2,:,:,izp)
    )
    vx = vx*hbar;
    vy = vy*hbar;
    vz = vz*hbar;
    vx,vy,vz
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
    psi.*exp(1.0im*alpha);
end

"""
multiplies exp(i*q) to (psi1,psi2)
"""
function GaugeTransform(psi1, psi2, q)
    eiq = exp(1.0im*q);
    psi1 = psi1.*eiq
    psi2 = psi2.*eiq
    psi1,psi2
end


"""
extracts Clebsch variable s=(sx,sy,sz) from (psi1,psi2)
"""
function Hopf(psi1,psi2)
    a = real(psi1)
    b = imag(psi1)
    c = real(psi2)
    d = imag(psi2)
    sx = 2*(a.*c + b.*d);
    sy = 2*(a.*d - b.*c);
    sz = a.^2 + b.^2 - c.^2 - d.^2;
    sx,sy,sz
end

"""
normalizes (psi1,psi2)
"""
function Normalize(psi1,psi2)
    psi_norm = sqrt(abs(psi1).^2 + abs(psi2).^2);
    psi1 = psi1./psi_norm;
    psi2 = psi2./psi_norm;
    psi1, psi2
end
