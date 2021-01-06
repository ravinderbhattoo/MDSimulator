export acceleration!, potential_energy

# acceleration function for MLPotentials
function MDBase.acceleration!(dv, v::Array{F,2}, u::Array{F,2}, pots::Array{M,1}, params) where {F <: Number, M <: MLModelPotential}
    for pot in pots
        dv .+= acceleration(v, u, pot)
    end
end


# energy function for MLPotentials
function MDBase.potential_energy(v::Array{F,2}, u::Array{F,2}, pots::Array{M,1}, params) where {F <: Number, M <: MLModelPotential}
    pe = 0.0
    for pot in pots
        pe += potential_energy(v, u, pot)
    end
    return pe
end
