export cutoff, reneighboring!

cutoff(pot::Any) = pot.R
cutoff(pot::MLPotential) = 0.0UNITS.distance

function reneighboring!(x, integrator, params)
    neigh_mat = find_neighbors(x, params.S.others.cut_dist, params.S.sim.boundary_condition, params.S.sim.others.max_neighs, hard_max = params.S.sim.others.max_neighs_hard_set)
    rows = size(neigh_mat, 1)
    fill!(params.S.others.neighs, 0)
    nrows = min(rows, size(params.S.others.neighs, 1))
    if nrows < rows
        print("Neighs has exceeded($rows) capacity($nrows)\n")
    end
    params.S.others.neighs[1:nrows,:] = neigh_mat[1:nrows,:]
end
