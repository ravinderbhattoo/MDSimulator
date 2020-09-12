export get_acceleration, get_potential_energy


function get_acceleration(v, u, params)
    dv = zeros(size(u))
    for pot in params.interatomic_potentials
        if !(typeof(pot) <: MLPotential)
            acceleration!(dv, u, pot, params.a_ids, params.m_ids, params.mass, params.boundary_condition, params.neighs)
        else
            acceleration!(dv, v, u, pot)
        end
    end
    dv *= params.uf
end


function get_potential_energy(u, params)
    e = 0.0
    for pot in params.interatomic_potentials
        if !(typeof(pot) <: MLPotential)
            e += potential_energy(u, pot, params.a_ids, params.m_ids, params.mass, params.boundary_condition, params.neighs)
        else
            nothing
        end
    end
    e
end

include("./pair.jl")
include("./ML.jl")
