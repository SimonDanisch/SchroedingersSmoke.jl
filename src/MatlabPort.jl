# Baseline version:
module MatlabPort

# matlab compat functions
_land(a, b) = a != 0 && b != 0
function land(a, b)
    map(_land, a, b)
end

ndgrid(v::AbstractVector) = copy(v)
function ndgrid(v1::AbstractVector, v2::AbstractVector)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end
function ndgrid{T}(vs::AbstractVector{T}...)
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array(T, sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

type ISF
    px::Array{Float64,3}; py::Array{Float64,3}; pz::Array{Float64,3}          # coordinates of grid points
    ix::UnitRange{Int}; iy::UnitRange{Int}; iz::UnitRange{Int}          # 1D index array
    iix::Array{Int, 3}; iiy::Array{Int, 3}; iiz::Array{Int, 3}         # 3D index array
    dx::Float64; dy::Float64; dz::Float64          # edge length
    sizex::Int; sizey::Int; sizez::Int # size of grid
    resx::Int; resy::Int; resz::Int    # number of grid points in each dimension

    hbar::Float64             # reduced Planck constant
    dt::Float64               # time step
    SchroedingerMask::Array{Complex{Float64},3} # Fourier coefficient for solving Schroedinger eq


    function ISF(vol_size::NTuple{3}, vol_res::NTuple{3}, hbar, dt)
        obj = new()
        obj.sizex, obj.sizey, obj.sizez = vol_size
        obj.resx = round(Int, vol_res[1])
        obj.resy = round(Int, vol_res[2])
        obj.resz = round(Int, vol_res[3])
        obj.dx = obj.sizex / obj.resx
        obj.dy = obj.sizey / obj.resy
        obj.dz = obj.sizez / obj.resz
        obj.ix = 1:obj.resx
        obj.iy = 1:obj.resy
        obj.iz = 1:obj.resz
        obj.iix, obj.iiy, obj.iiz = ndgrid(obj.ix,obj.iy,obj.iz)
        obj.px = (obj.iix-1) * obj.dx
        obj.py = (obj.iiy-1) * obj.dy
        obj.pz = (obj.iiz-1) * obj.dz
        obj.hbar = hbar
        obj.dt = dt
        BuildSchroedinger(obj)
        obj
    end

end

function BuildSchroedinger(obj::ISF)
    nx=obj.resx; ny=obj.resy; nz=obj.resz;
    fac = -4*pi^2*obj.hbar;
    kx = (obj.iix-1-nx/2) / (obj.sizex);
    ky = (obj.iiy-1-ny/2) / (obj.sizey);
    kz = (obj.iiz-1-nz/2) / (obj.sizez);
    lambda = fac*(kx.^2+ky.^2+kz.^2);
    obj.SchroedingerMask = exp(1.0im*lambda*obj.dt/2.);
end



function VelocityOneForm(obj, psi1, psi2, hbar=1.0)
    ixp = mod(obj.ix, obj.resx) + 1;
    iyp = mod(obj.iy, obj.resy) + 1;
    izp = mod(obj.iz, obj.resz) + 1;
    vx = angle(
        conj(psi1).*view(psi1, ixp,:,:) +
        conj(psi2).*view(psi2, ixp,:,:)
    );
    vy = angle(
        conj(psi1).*view(psi1,:,iyp,:) +
        conj(psi2).*view(psi2,:,iyp,:)
    )
    vz = angle(
        conj(psi1).*view(psi1,:,:,izp) +
        conj(psi2).*view(psi2,:,:,izp)
    )
    vx = vx*hbar;
    vy = vy*hbar;
    vz = vz*hbar;
    vx,vy,vz
end

function PoissonSolve(obj, f)
    f = fft(f)
    sx = sin(pi*(obj.iix-1)/obj.resx)/obj.dx
    sy = sin(pi*(obj.iiy-1)/obj.resy)/obj.dy
    sz = sin(pi*(obj.iiz-1)/obj.resz)/obj.dz
    denom = sx.^2 + sy.^2 + sz.^2
    fac = -0.25./denom
    fac[1,1,1] = 0.0
    f = f .* fac
    ifft(f)
end

function GaugeTransform(psi1, psi2, q)
    eiq = exp(1.0im * q);
    psi1 = psi1 .* eiq
    psi2 = psi2 .* eiq
    psi1, psi2
end

function PressureProject(obj, psi1, psi2)
    vx,vy,vz = VelocityOneForm(obj, psi1,psi2)
    div = Div(obj, vx,vy,vz)
    q = -PoissonSolve(obj, div)
    x = GaugeTransform(psi1, psi2, q)
    x
end

function Normalize(psi1, psi2)
    psi_norm = sqrt(abs(psi1).^2 + abs(psi2).^2);
    psi1 = psi1 ./ psi_norm
    psi2 = psi2 ./ psi_norm
    psi1, psi2
end

function Div(obj, vx, vy, vz)
    ixm = mod(obj.ix-2, obj.resx) + 1
    iym = mod(obj.iy-2, obj.resy) + 1
    izm = mod(obj.iz-2, obj.resz) + 1
    f =     (vx - view(vx, ixm,:,:))/(obj.dx^2)
    f = f + (vy - view(vy, :,iym,:))/(obj.dy^2)
    f = f + (vz - view(vz, :,:,izm))/(obj.dz^2)
    f
end

end
