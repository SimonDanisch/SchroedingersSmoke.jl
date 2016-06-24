
type TorusDEC
    p::cl.CLArray{Point3f0,3}          # coordinates of grid points
    fac::cl.CLArray{Float32, 3}
    ix::UnitRange{Int}; iy::UnitRange{Int}; iz::UnitRange{Int}          # 1D index array
    iix::Array{Int, 3}; iiy::Array{Int, 3}; iiz::Array{Int, 3}         # 3D index array
    d::Vector{Float32}          # edge length
    size::Vector{Int32} # size of grid
    res::Vector{Int32}    # number of grid points in each dimension
    velocity::Vector{Vec3f0}

    function TorusDEC(vol_size::NTuple{3}, vol_res::NTuple{3})
        obj = new()
        obj.sizex,obj.sizey, obj.sizez = vol_size
        obj.resx = round(Int, vol_res[1])
        obj.resy = round(Int, vol_res[2])
        obj.resz = round(Int, vol_res[3])
        obj.dx = obj.sizex/obj.resx
        obj.dy = obj.sizey/obj.resy
        obj.dz = obj.sizez/obj.resz

        ix, iy, iz = 1:obj.resx, 1:obj.resy, 1:obj.resz

        obj.iix, obj.iiy, obj.iiz = ndgrid(ix, iy, iz)
        obj.px = cl.CLArray((obj.iix-1)*obj.dx)
        obj.py = cl.CLArray((obj.iiy-1)*obj.dy)
        obj.pz = cl.CLArray((obj.iiz-1)*obj.dz)

        sx = sin(pi*(obj.iix-1)/obj.resx)/obj.dx
        sy = sin(pi*(obj.iiy-1)/obj.resy)/obj.dy
        sz = sin(pi*(obj.iiz-1)/obj.resz)/obj.dz
        denom = sx.^2 + sy.^2 + sz.^2
        fac = -0.25./denom
        fac[1,1,1] = 0.0
        obj.fac = cl.CLArray(fac)
        obj.velocity = zeros(Vec3f0, vol_size)
        
        return obj
    end

end

@cl_kernel Div program_torusDEC velocity f res d


@cl_kernel StaggeredSharp program_torusDEC velocity res d_inv

function StaggeredSharp(obj::TorusDEC)
    StaggeredSharp(obj.velocity, StaggeredSharp.res, 1f0/obj.d)
end

"""
PoissonSolve by Spectral method
"""
function PoissonSolve(obj, f)
    f = fft(f)
    mul!(f, obj.fac)
    ifft(f)
end
