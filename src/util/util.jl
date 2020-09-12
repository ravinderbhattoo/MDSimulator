export cutoff, reneighboring!

cutoff(pot::Any) = pot.R
cutoff(pot::MLPotential) = 0.0UNITS.distance

function reneighboring!(x, integrator, params)
    neigh_mat = find_neighbors(x, params.cut_dist, params.boundary_condition, params.max_neighs, hard_max = params.max_neighs_hard_set)
    rows = size(neigh_mat, 1)
    fill!(integrator.p[2].neighs, 0)
    nrows = min(rows, size(integrator.p[2].neighs)[1])
    if nrows < rows
        print("Neighs has exceeded($rows) capacity($nrows)\n")
    end
    integrator.p[2].neighs[1:nrows,:] = neigh_mat[1:nrows,:]
end
