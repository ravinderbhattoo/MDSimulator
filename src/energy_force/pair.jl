export acceleration!, potential_energy

# acceleration function for PairPotentials
function MDBase.acceleration!(dv, v, u, pot::T, params) where T <: PairPotential
    neighs = params.S.others.neighs
    BC = params.S.sim.boundary_condition
    for i in 1:size(u, 2)
        for j in 2:neighs[1,i]+1
            r, r2, dr = distance(u[:,i], u[:,neighs[j,i]], BC)
            mag = potential_force(r, pot)/(params.S.sim.mass[i]*r)
            dv[1,i] += dr[1]*mag
            dv[2,i] += dr[2]*mag
            dv[3,i] += dr[3]*mag
        end
    end
end


# energy function for PairPotentials
function MDBase.potential_energy(u, pot::T, params) where T <: PairPotential
    e = 0
    neighs = params.S.others.neighs
    for i in 1:size(u, 2)
        for j in 2:neighs[1,i]+1
            r, r2, dr = distance(u[:,i], u[:,neighs[j,i]], params.S.sim.boundary_condition)
            mag = potential_energy(r, pot)
            e += mag
        end
    end
    return e/2
end
