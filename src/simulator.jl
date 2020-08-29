using DifferentialEquations
using ParticlesTools
using Unitful

# exports
export MDSim, simulate

struct MDSim
    units
    initial_condition
    mass
    interatomic_potentials
    boundary_condition
    a_ids
    m_ids
    Δτ
    save_every
    global_save_every
    max_neighs
    max_neighs_hard_set
    reneighboring_every
end

function MDSim(initial_condition, mass, interatomic_potentials, boundary_condition; a_ids = nothing, m_ids = nothing, Δτ = 1.0e-3, save_every = 1000, global_save_every = 100, max_neighs = 200, max_neighs_hard_set = 100, reneighboring_every = 100)

    if a_ids==nothing
        a_ids = ones(size(initial_condition)[1][2])
        m_ids = copy(a_ids)
    elseif m_ids==nothing
        m_ids = copy(a_ids)
    end

    # Units
    u_x, u_v = [unit(i[1]) for i in initial_condition]
    units = (u_x=u_x, u_v=u_v, u_mass=unit(mass), u_f_interatomic_potentials=[unit(potential_force(1u_x)) for pot in interatomic_potentials])

    MDSim(units, initial_condition, mass, interatomic_potentials, boundary_condition, a_ids, m_ids, Δτ, save_every, global_save_every, max_neighs, max_neighs_hard_set, reneighboring_every)
end


function MD_soode(v, u, p, t)
    emvironment, params = p
    dv = get_acceleration(u, params)
    # for env in environment
    #     @. dv += ddu(env)
    # end
    return dv
end

function simulate(Sim::MDSim, n::Int64, environment)
    tspan = (0.0, n*Sim.Δτ)

    u0, v0 = Sim.initial_condition

    cut_dist = maximum([pot.R for pot in Sim.interatomic_potentials])

    neighs = find_neighbors(u0, cut_dist, Sim.boundary_condition, Sim.max_neighs, hard_max=Sim.max_neighs_hard_set)

    nn = Int64(n/Sim.global_save_every)
    global_vals = [WellArray(nn), WellArray(nn), WellArray(nn)]


    params = (interatomic_potentials=Sim.interatomic_potentials, a_ids=Sim.a_ids, m_ids=Sim.m_ids, mass=Sim.mass, boundary_condition=Sim.boundary_condition, neighs=neighs, uf=Sim.unit_factor, N=size(u0)[2], cut_dist=cut_dist, max_neighs=Sim.max_neighs, max_neighs_hard_set=Sim.max_neighs_hard_set, reneighboring_every=Sim.reneighboring_every, global_save_every=Sim.global_save_every, global_vals=global_vals)

    p = (environment, params)

    soode(v, u, p, t) = MD_soode(v, u, p, t)

    prob = SecondOrderODEProblem(soode, v0, u0, tspan, p, callback=DiscreteCallback(condition_f, affect_f!, save_positions=(false,false)))

    solve(prob, VelocityVerlet(), dt=Sim.Δτ, saveat=0.0:Sim.save_every*Sim.Δτ:n*Sim.Δτ)
end


function get_acceleration(u, params)
    dv = zeros(size(u))
    for pot in params.interatomic_potentials
        dv += params.uf*acceleration(u, pot, params.a_ids, params.m_ids, params.mass, params.boundary_condition, params.neighs)
    end
    dv
end


function get_potential_energy(u, params)
    e = 0.0
    for pot in params.interatomic_potentials
        e += potential_energy(u, pot, params.a_ids, params.m_ids, params.mass, params.boundary_condition, params.neighs)
    end
    e
end


function condition_f(u,t,integrator)
    environment, params = integrator.p
    N = params.N
    v = reshape(Array(integrator.u),(3,:))[:,1:N]
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    apply_simulation_bc!(x, v, params.boundary_condition)
    integrator.u[:] = hcat(reshape(v,(1,:)),reshape(x,(1,:)))
    integrator.sol.u[end][:] = integrator.u[:]
    if round(integrator.t/integrator.dt)%params.global_save_every==0
        fillit!(params.global_vals[1], 0.5*sum(params.mass .* sum(v.*v, dims=1)))
        fillit!(params.global_vals[3], get_potential_energy(x, params))
        print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \tKE: $(params.global_vals[1][end]) \tPE: $(params.global_vals[3][end]) \n")
    end
    isapprox(round(integrator.t/integrator.dt)%params.reneighboring_every,0.0)
end

function affect_f!(integrator)
    params = integrator.p[2]
    N = params.N
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \nReneighboring(cutoff = $(params.cut_dist))...")
    neigh_mat = find_neighbors(x, params.cut_dist, params.boundary_condition, params.max_neighs, hard_max = params.max_neighs_hard_set)
    rows = size(neigh_mat)[1]
    fill!(integrator.p[2].neighs, 0)
    nrows = min(rows, size(integrator.p[2].neighs)[1])
    if nrows < rows
        print("Neighs has exceeded($rows) capacity($nrows)\n")
    end
    integrator.p[2].neighs[1:nrows,:] = neigh_mat[1:nrows,:]
    print("\tDone\n")
end
