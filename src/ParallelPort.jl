module ParallelPort

using GeometryTypes
using StaticArrays
using GPUArrays
const ArrayType = GPUArrays.JLArray

# lots of parameters. To lazy to write out types,
# still don't want to waste performance
type ISF{IntType, FloatType}

    grid_res::Vec{3, IntType}
    physical_size::Vec{3, IntType}
    d::Vec{3, FloatType}
    velocity::ArrayType{Vec{3, FloatType}, 3}
    # gridpoints in pysical size
    ranges::Vec{3, Vector{FloatType}}

    hbar::FloatType             # reduced Planck constant
    dt::FloatType               # time step
    mask::ArrayType{Complex{FloatType}, 3} # Fourier coefficient for solving Schroedinger eq
    # precalculated values
    idx_shifted::Vec{3, Vector{IntType}} # indices shifted by one
    fac::ArrayType{Complex{FloatType}, 3}
    # temporary allocations for intermediates
    f::ArrayType{FloatType, 3}

    function ISF(physical_size, dims, hbar, dt)

        VT = Vec{3, FloatType}
        grid_res = Vec(dims)
        physical_size = Vec(physical_size)
        d = VT(physical_size ./ grid_res)
        fac = zeros(FloatType, dims)
        mask = zeros(Complex{FloatType}, dims)
        f = -4 * pi^2 * hbar
        for x = 1:dims[1], y = 1:dims[2], z = 1:dims[3]
            xyz = Vec(x, y, z)
            # fac
            s = sin(pi * (xyz - 1) ./ grid_res) ./ d
            denom = sum(s .^ 2)
            fac[x,y,z] = -0.25f0 ./ denom

            # schroedingers mask
            k = (xyz - 1 - grid_res ./ 2) ./ physical_size
            lambda = f * sum(k .^ 2)
            mask[x,y,z] = exp(1.0im * lambda * dt / 2)
        end
        fac[1,1,1] = 0

        f = zeros(FloatType, dims)
        velocity = zeros(VT, dims)

        ranges = Vec(ntuple(3) do i
             ((1:dims[i]) - 1) * d[i]
        end)
        idx_shifted = Vec(ntuple(3) do i
            mod(1:grid_res[i], grid_res[i]) + 1
        end)
        new(
            grid_res,
            physical_size,
            d,
            velocity,
            ranges,
            hbar,
            dt,
            mask,
            idx_shifted,
            fac,
            f
        )
    end
end


import BroadcastPort: velocity_one_form!
import BroadcastPort: gauge_transform!
import BroadcastPort: poisson_solve!
import BroadcastPort: pressure_project!
import BroadcastPort: Normalize!
import BroadcastPort: schroedinger_flow!
import BroadcastPort: Particles!
import BroadcastPort: staggered_advect!


end
