using GPUArrays, CUDAnative
import GPUArrays: GPUArray, CUBackend, cu_map
cuctx = CUBackend.init()
using GeometryTypes

function VelocityOneForm(obj, psi, hbar=1.0f0)
    map_idx!(obj.t.velocity, psi) do idx, psi
        x,y,z = idx.I
        ixp, iyp, izp = mod(Vec{3, Int}(idx.I), obj.t.res) + 1
        @inbounds begin
            psi1,  psi2  = psi[x,y,z]
            psix1, psix2 = psi[ixp,y,z]
            psiy1, psiy2 = psi[x,iyp,z]
            psiz1, psiz2 = psi[x,y,izp]
        end
        psi1c = conj(psi1); psi2c = conj(psi2);
        vx = angle(psi1c*psix1 + psi2c*psix2)
        vy = angle(psi1c*psiy1 + psi2c.*psiy2)
        vz = angle(psi1c*psiz1 + psi2c*psiz2)
        Vec3f0(vx, vy, vz)*hbar
    end
end
function GaugeTransform(psi, q)
    broadcast!(psi, psi, q) do psi, q
        eiq = exp(1f0*im*q)
        (psi[1]*eiq, psi[2]*eiq)
    end
    psi
end
function Div(obj::TorusDEC, velocity, f)
    map_idx!(f, velocity) do idx, velocity
        x,y,z = idx.I # cartesian index
        ixm = mod(x-2, obj.resx) + 1
        iym = mod(y-2, obj.resy) + 1
        izm = mod(z-2, obj.resz) + 1
        _x = velocity[ixm, y, z][1]
        _y = velocity[x, iym, z][2]
        _z = velocity[x, y, izm][3]
        v  = velocity[x,y,z]
        ff =  (v[1] - _x) / obj.dx^2
        ff += (v[2] - _y) / obj.dy^2
        ff += (v[3] - _z) / obj.dz^2
        ff
    end
    f
end
function PoissonSolve(obj, f)
    fc = fft(f)
    fc .= (*).(fc, obj.fac)
    ifft!(fc)
    fc
end
function PressureProject(obj, psi)
    velocity = VelocityOneForm(obj, psi)
    div = Div(obj.t, velocity)
    q = PoissonSolve(obj.t, div);
    GaugeTransform(psi1, psi2, (-).(q))
end

n = 50
psi_1_2 = GPUArray([(rand(Complex64), rand(Complex64)) for i=1:n, j=1:n, k=1:n]);
s = GPUArray(zeros(Vec3f0, (n,n,n)))
Hopf(psi_1_2, s)
