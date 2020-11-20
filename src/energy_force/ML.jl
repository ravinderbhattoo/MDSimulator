export acceleration!, potential_energy

# acceleration function for MLPotentials
function MDBase.acceleration!(dv, v::Array{F,2}, u::Array{F,2}, pot::M, params) where {F <: Number, M <: MLModelPotential}
    dv .+= acceleration(v, u, pot)
end


# energy function for MLPotentials
function MDBase.potential_energy(v::Array{F,2}, u::Array{F,2}, pot::M, params) where {F <: Number, M <: MLModelPotential}
    return potential_energy(v, u, pot)
end
