tuple_dot(a, b) = sum(a .+ b)

"""
function returning true at nozzle position
"""
function isjet(p, nozzle_cen, nozzle_len, nozzle_rad)
    (abs(p[1] - nozzle_cen[1]) <= nozzle_len / 2) &
    ((p[2] - nozzle_cen[2])^2 +
    (p[3] - nozzle_cen[3])^2 <= nozzle_rad .^ 2) != 0
end

function restrict_kernel(psi, isjet, kvec, pos, omgterm)
    if isjet
        amp = @fastmath abs.(psi)
        phase = tuple_dot(kvec, pos) - omgterm
        @fastmath amp .* exp(Complex64(0f0, 1f0) .* phase)
    else
        psi
    end
end

"""
Restrict the velocity in restrict_positions
"""
function restrict_velocity!(isf, psi, kvec, restrict_positions, omgterm = 0f0)
    psi .= restrict_kernel.(psi, restrict_positions, (kvec,), isf.positions, omgterm)
end



"""
Add particles around a jet on the CPU
"""
function add_particles!(particles::Vector, range, nozzle_cen, nozzle_rad)
    randstate = GPUArrays.cached_state(particles)
    for i in range
        rt = rand(Float32) * 2f0 * pi

        spawn_offset = ((
            rand(Float32),
            rand(Float32),
            rand(Float32),
        ) .* 0.1f0) .- 0.05f0

        spawned = ((
            nozzle_cen[1] - 0.1f0,
            nozzle_cen[2] + 0.9f0 * nozzle_rad * cos(rt),
            nozzle_cen[3] + 0.9f0 * nozzle_rad * sin(rt)
        )) .+ spawn_offset

        particles[i] = spawned
    end
    return
end


"""
Add particles around a jet on the GPU!
"""
function add_particles!(particles, range, nozzle_cen, nozzle_rad)
    randstate = GPUArrays.cached_state(particles)
    N = length(range)
    gpu_call(
        add_particles!, particles,
        (
            randstate, particles, UInt32(first(range)), UInt32(N),
            nozzle_cen, nozzle_rad, Float32(pi)
        ),
        N
    )
    return
end


function add_particles!(state, randstate, particles, offset, n, nozzle_cen, nozzle_rad, pi)
    ilin = linear_index(state)
    ilin > n && return
    rt = gpu_rand(Float32, state, randstate) * 2f0 * pi

    spawn_offset = ((
        gpu_rand(Float32, state, randstate),
        gpu_rand(Float32, state, randstate),
        gpu_rand(Float32, state, randstate),
    ) .* 0.1f0) .- 0.05f0

    spawned = ((
        nozzle_cen[1] - 0.1f0,
        nozzle_cen[2] + 0.9f0 * nozzle_rad * cos(rt),
        nozzle_cen[3] + 0.9f0 * nozzle_rad * sin(rt)
    )) .+ spawn_offset

    particles[ilin + offset] = spawned
    return
end


export add_particles!, restrict_velocity!, isjet
