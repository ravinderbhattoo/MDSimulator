export acceleration!, potential_energy!

# acceleration function for PairPotentials
function acceleration!(a, x, pot::PairPotential, a_ids, m_ids, mass, BC, neighs)
    for i in 1:size(x, 2)
        for j in 2:neighs[1,i]+1
            r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
            mag = potential_force(r, pot)/(mass[i]*r)
            a[1,i] += dr[1]*mag
            a[2,i] += dr[2]*mag
            a[3,i] += dr[3]*mag
        end
    end
end


# energy function for PairPotentials
function potential_energy(x, pot::PairPotential, a_ids, m_ids, mass, BC, neighs)
    e = 0
    for i in 1:size(x, 2)
        for j in 2:neighs[1,i]+1
            r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
            mag = potential_energy(r, pot)
            e += mag
        end
    end
    return e/2
end

# energy function (per atom) for PairPotentials
function potential_energy_per_atom!(pe, x, pot::PairPotential, a_ids, m_ids, mass, BC, neighs)
    fill!(pe, 0.0)
    for i in 1:size(x, 2)
        for j in 2:neighs[1,i]+1
            r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
            mag = potential_energy(r, pot)
            pe[i] += mag
        end
    end
end
