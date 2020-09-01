
export acceleration, potential_energy

# acceleration function for PairPotentials
function acceleration(x, pot::PairPotential, a_ids, m_ids, mass, BC, neighs)
    a = zeros(size(x))
    for i in 1:size(x)[2]
        for j in 1:size(neighs)[1]
            if neighs[j,i]==0
                break
            else
                r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
                mag = potential_force(r, pot)/(mass[i]*r)
                a[1,i] += dr[1]*mag
                a[2,i] += dr[2]*mag
                a[3,i] += dr[3]*mag
            end
        end
    end
    return a
end


# energy function for PairPotentials
function potential_energy(x, pot::PairPotential, a_ids, m_ids, mass, BC, neighs)
    e = 0
    for i in 1:size(x, 2)
        for j in 1:size(neighs, 1)
            if neighs[j,i]==0
                break
            else
                r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
                mag = potential_energy(r, pot)
                e += mag
            end
        end
    end
    return e/2
end
