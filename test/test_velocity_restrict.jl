using GeometryTypes

# straight matlab port
function VelocityOneForm(psi1, psi2, hbar=1.0)
    dims = size(psi1)
    ix, iy, iz = 1:dims[1], 1:dims[2], 1:dims[3]
    ixp = mod(ix, dims[1]) + 1
    iyp = mod(iy, dims[2]) + 1
    izp = mod(iz, dims[3]) + 1
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
    map(Vec3f0, zip(vx,vy,vz))
end

function map_idx!{F, T}(f::F, A::AbstractArray, args::T)
    n = length(A); dims = size(A)
    Base.Threads.@threads for i = 1:n
        idx = Vec{3, Int}(ind2sub(dims, i))
        val = f(idx, A, args)
        @inbounds A[i] = val
    end
    A
end
@inline function inner_velocity_one_form(i, velocity, idx_psi_hbar)
    idx2, psi, hbar = idx_psi_hbar
    i2 = (idx2[1][i[1]], idx2[2][i[2]], idx2[3][i[3]])
    @inbounds begin
        psi12  = psi[i[1],  i[2],  i[3]]
        psix12 = psi[i2[1], i[2],  i[3]]
        psiy12 = psi[i[1],  i2[2] ,i[3]]
        psiz12 = psi[i[1],  i[2],  i2[3]]
    end
    psi1n = Vec(psix12[1], psiy12[1], psiz12[1])
    psi2n = Vec(psix12[2], psiy12[2], psiz12[2])
    angle.(
        conj(psi12[1]) .* psi1n .+
        conj(psi12[2]) .* psi2n
    ) * hbar
end
function velocity_one_form!(velocity, psi)
    dims = size(psi)
    idx = ntuple(3) do i
        mod(1:dims[i], dims[i]) + 1
    end
    arg = (idx, psi, 1f0)
    map_idx!(inner_velocity_one_form, velocity, arg)
end

dims = (64, 32, 32)
psi1, psi2 = rand(Complex64, dims), rand(Complex64, dims);
psi = map(identity, zip(psi1, psi2));
velocity = zeros(Vec3f0, dims);

v1 = velocity_one_form!(velocity, psi);
v2 = VelocityOneForm(psi1, psi2);

all(isapprox.(v1, v2))
