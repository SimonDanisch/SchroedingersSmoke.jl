type TorusDEC
    p::CLArray{Point3f0,3}          # coordinates of grid points
    fac::CLArray{Float32, 3}
    ix::UnitRange{Int}; iy::UnitRange{Int}; iz::UnitRange{Int}          # 1D index array
    iix::AFArray{Int32, 3}; iiy::AFArray{Int32, 3}; iiz::AFArray{Int32, 3}         # 3D index array
    d::Vector{Float32}          # edge length
    size::Vector{Int32} # size of grid
    res::Vector{Int32}    # number of grid points in each dimension

    function TorusDEC(vol_size::NTuple{3}, vol_res::NTuple{3})
        obj = new()
        obj.sizex,obj.sizey, obj.sizez = vol_size
        obj.resx = round(Int, vol_res[1])
        obj.resy = round(Int, vol_res[2])
        obj.resz = round(Int, vol_res[3])
        obj.dx = obj.sizex/obj.resx
        obj.dy = obj.sizey/obj.resy
        obj.dz = obj.sizez/obj.resz
        obj.ix = 1:obj.resx
        obj.iy = 1:obj.resy
        obj.iz = 1:obj.resz
        obj.iix,obj.iiy,obj.iiz = ndgrid(obj.ix,obj.iy,obj.iz)
        obj.px = CLArray((obj.iix-1)*obj.dx)
        obj.py = CLArray((obj.iiy-1)*obj.dy)
        obj.pz = CLArray((obj.iiz-1)*obj.dz)

        sx = sin(pi*(obj.iix-1)/obj.resx)/obj.dx
        sy = sin(pi*(obj.iiy-1)/obj.resy)/obj.dy
        sz = sin(pi*(obj.iiz-1)/obj.resz)/obj.dz
        denom = sx.^2 + sy.^2 + sz.^2
        fac = -0.25./denom
        fac[1,1,1] = 0.0
        obj.fac = CLArray(fac)
        return obj
    end
    function TorusDEC(varargin...)
        n = length(varargin)
        if n == 0 # empty instance
            return
        elseif n == 4
            mi = findmax(varargin[1:3])
            ratio = [varargin[1:3]...]./varargin[mi]
            res = round(ratio*varargin[4])
            return TorusDEC(varargin[1:3],res[1],res[2],res[3])
        end
        error(
            "TorusDEC:badinput
            Wrong number of inputs."
        )
    end
end



"""
For a 1-form v compute the function `*d*v`
"""
@cl_kernel Div velocity res d



"""
 For a 1-form v compute the corresponding vector field `v^sharp` as
 a staggered vector field living on edges
"""
@cl_kernel StaggeredSharp velocity d_inv res
function StaggeredSharp(obj::TorusDEC)
    StaggeredSharp(obj.velocity, StaggeredSharp.res, 1f0/obj.d)
end

"""
PoissonSolve by Spectral method
"""
function PoissonSolve(obj, f)
    f = fft(f)
    f = f .* fac
    ifft(f)
end
