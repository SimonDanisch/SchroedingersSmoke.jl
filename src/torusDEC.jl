
type TorusDEC
    p::cl.CLArray{Point3f0,3}          # coordinates of grid points
    fac::cl.CLArray{Float32, 3}
    iix::Array{Int, 3}; iiy::Array{Int, 3}; iiz::Array{Int, 3}         # 3D index array
    d::Vector{Float32}          # edge length
    size::Vector{Int32} # size of grid
    res::Vector{Int32}    # number of grid points in each dimension
    velocity::Vector{Vec3f0}

    function TorusDEC(vol_size::NTuple{3}, vol_res::NTuple{3})
        obj = new()
        obj.size = Int32[vol_size...]
        obj.res = Int32[round(Int, vol_res[i]) for i=1:3]
        obj.d = Float32[obj.size[i]/obj.res[i] for i=1:3]

        ix, iy, iz = ntuple(i->1:obj.res[i], 3)

        obj.iix, obj.iiy, obj.iiz = ndgrid(ix, iy, iz)
        obj.px = cl.CLArray((obj.iix-1)*obj.d[1])
        obj.py = cl.CLArray((obj.iiy-1)*obj.d[2])
        obj.pz = cl.CLArray((obj.iiz-1)*obj.d[3])

        sx = sin(pi*(obj.iix-1)/obj.resx)/obj.d[1]
        sy = sin(pi*(obj.iiy-1)/obj.resy)/obj.d[2]
        sz = sin(pi*(obj.iiz-1)/obj.resz)/obj.d[3]
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
