export summary_, get_kinetic_energy, maxwell_boltzmann_velocity, get_temperature, get_momentum_all

function get_temperature(v::Array{T1,2}, m::Array{T2, 1}) where {T1, T2}
    ke = get_kinetic_energy(v, m)
    N = length(m)
    return 2ke/(3N*CONSTANTS.kb)
end

function get_temperature(ke::Number, N::Int)
    return 2ke/(N*CONSTANTS.kb)
end

function maxwell_boltzmann_velocity(N, T, m)
    v_dev = uconvert(UNITS.velocity, sqrt(CONSTANTS.kb * T / m))
    v0x = v_dev*randn((1,N))
    v0y = v_dev*randn((1,N))
    v0z = v_dev*randn((1,N))
    v0 = [v0x; v0y; v0z]
    v0 .-= sum(v0, dims=2)/N
end

function get_kinetic_energy(v::Array{T1, 2}, m::Array{T2, 1}) where {T1, T2}
    return 0.5*sum((v.^2)*m)
end

function summary_(params)
    KE = params.S.others.thermo_vals.ke
    PE = params.S.others.thermo_vals.pe
    temp = params.S.others.thermo_vals.temp
    return KE, temp, PE
end


function get_momentum(res, t, params)
    res(t).x[1]*params.S.sim.mass
end

function get_momentum_all(res, params)
    p = zeros(3,length(res.t))
    for (ind, t) in enumerate(res.t)
        p[:,ind] = get_momentum(res, t, params)
    end
    return p
end

#
