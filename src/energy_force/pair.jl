export acceleration!, potential_energy

# acceleration function for PairPotentials
function MDBase.acceleration!(dv, v, u, pot::T, params) where T <: Union{PairPotential, MLPairPotential}
    neighs = params.S.others.neighs
    BC = params.S.sim.boundary_condition
    a_ids = params.S.sim.a_ids
    @inbounds for i in 1:size(u, 2)
        ax, ay, az = u[1,i], u[2,i], u[3,i]
        N = neighs[1,i]+1
        for j in 2:N
            k = neighs[j,i]
            bx, by, bz = u[1,k], u[2,k], u[3,k]
            r, r2, dr = distance((ax, ay,az), (bx, by, bz), BC)
            mag = potential_force(r, pot, (a_ids[i],a_ids[k]))/(params.S.sim.mass[i]*r)
            dv[1,i] += dr[1]*mag
            dv[2,i] += dr[2]*mag
            dv[3,i] += dr[3]*mag
        end
    end
end


# energy function for PairPotentials
function MDBase.potential_energy(v, u, pot::T, params) where T <: Union{PairPotential, MLPairPotential}
    e = 0.0
    neighs = params.S.others.neighs
    a_ids = params.S.sim.a_ids
    @inbounds for i in 1:size(u, 2)
        ax, ay, az = u[1,i], u[2,i], u[3,i]
        N = neighs[1,i]+1
        for j in 2:N
            k = neighs[j,i]
            bx, by, bz = u[1,k], u[2,k], u[3,k]
            r, r2, dr = distance((ax, ay,az), (bx, by, bz), params.S.sim.boundary_condition)
            mag = potential_energy(r, pot, (a_ids[i], a_ids[k]))
            e += mag
        end
    end
    return e/2
end
