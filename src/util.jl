export summary_, get_kinetic_energy, write_trajectory, save_trajectory, maxwell_boltzmann_velocity, get_temperature

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
    return 0.5*sum(m.*sum(v.^2, dims=1)')
end

function summary_(res)
    KE, temp, PE = res.prob.p[2].global_vals
    return KE, temp, PE
end


function write_trajectory(filename, res::DiffEqBase.ODESolution)
    if typeof(res) <: DiffEqBase.ODESolution
        for i in 1:length(res.t)
            t = res.t[i]
            data = Array(reshape(Array(res(t)),(3,:)))'
            N = Int(size(data, 1)/2)
            v = data[1:N,:]
            x = data[N+1:end,:]
            writedata("$(filename)_$i.data", (x, v))
        end
    else
        throw("Object is not of type DiffEqBase.ODESolution.")
    end
end


using JLD, HDF5
function save_trajectory(filename, res)
    jldopen("$filename", "w") do file
        g = g_create(file, "res") # create a group
        g["u"] = res.u              # create a scalar dataset inside the group
        g["t"] = res.t
        g["interatomic_potentials"] = res.prob.p[2].interatomic_potentials
        g["boundary_condition"] = res.prob.p[2].boundary_condition
        g["units"] = string(UNITS)
        g["global_vals"] = ustrip(res.prob.p[2].global_vals)
    end
    return nothing
end


#
