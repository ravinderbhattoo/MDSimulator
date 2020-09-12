export acceleration!, potential_energy!

# acceleration function for MLPotentials
function acceleration!(a, v, u, pot::MLPotential)
    a[:,:] += acceleration(v, u, pot)
end


# energy function for MLPotentials
function potential_energy(x, pot::MLPotential)
    return 0.0
end
