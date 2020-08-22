# exports
using DifferentialEquations
using ParticlesTools
export MDSim, simulate

struct MDSim
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
    MDSim(initial_condition, mass, interatomic_potentials, boundary_condition, a_ids, m_ids, Δτ, save_every, global_save_every, max_neighs, max_neighs_hard_set, reneighboring_every)
end


function MD_soode(v, u, p, t)
    emvironment, ode_params, cb_params = p
    dv = get_acceleration(u, ode_params)
    # for env in environment
    #     @. dv += ddu(env)
    # end
    return dv
end

function simulate(Sim::MDSim, n::Int64, environment)
    tspan = (0.0, n*Sim.Δτ)

    u0, v0 = Sim.initial_condition

    cut_dist = maximum([pot.R for pot in Sim.interatomic_potentials])

    neighs = find_neighbors(u0, cut_dist, Sim.max_neighs, hard_max=Sim.max_neighs_hard_set)

    ode_params = (Sim.interatomic_potentials, Sim.a_ids, Sim.m_ids, Sim.mass, Sim.boundary_condition, neighs)
    nn = Int64(n/Sim.global_save_every)
    global_vals = [WellArray(nn), WellArray(nn), WellArray(nn)]

    cb_params = (size(u0)[2], cut_dist, Sim.max_neighs, Sim.max_neighs_hard_set, Sim.reneighboring_every, Sim.global_save_every, global_vals)

    p = (environment, ode_params, cb_params)

    soode(v, u, p, t) = MD_soode(v, u, p, t)

    prob = SecondOrderODEProblem(soode, v0, u0, tspan, p, callback=DiscreteCallback(condition_f, affect_f!, save_positions=(false,false)))

    solve(prob, VelocityVerlet(), dt=Sim.Δτ, saveat=0.0:Sim.save_every*Sim.Δτ:n*Sim.Δτ)
end


function get_acceleration(u, params)
    pots, a_ids, m_ids, mass, bc, neighs = params
    dv = zeros(size(u))
    for pot in pots
        dv += acceleration(u, pot, a_ids, m_ids, mass, bc, neighs)
    end
    dv
end

function acceleration(x, pot, a_ids, m_ids, mass, BC, neighs)
    a = zeros(size(x))
    for i in 1:size(x)[2]
        for j in 1:size(neighs)[1]
            if neighs[j,i]!=0
                r, r2, dr = distance(x[:,i], x[:,neighs[j,i]], BC)
                mag = potential_force(r, pot)/(mass[i]*r) # , a_ids[i], a_ids[neighs[j,i]], m_ids[i], m_ids[neighs[j,i]]
                a[1,i] += dr[1]*mag
                a[2,i] += dr[2]*mag
                a[3,i] += dr[3]*mag
            end
        end
    end
    return a
end

function condition_f(u,t,integrator)
    environment, ode_params, cb_params = integrator.p
    N = cb_params[1]
    v = reshape(Array(integrator.u),(3,:))[:,1:N]
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    apply_simulation_bc!(x, v, ode_params[end-1])
    integrator.u[:] = hcat(reshape(v,(1,:)),reshape(x,(1,:)))
    integrator.sol.u[end][:] = integrator.u[:]
    if round(integrator.t/integrator.dt)%cb_params[6]==0
        fillit!(cb_params[7][1], 0.5*sum(ode_params[4] .* sum(v.*v, dims=1)))
        print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \tKE: $(cb_params[7][1][end])\n")
    end
    isapprox(round(integrator.t/integrator.dt)%cb_params[5],0.0)
end

function affect_f!(integrator)
    N, cut_dist, max_neighs, max_neighs_hard_set, _ = integrator.p[end]
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \nReneighboring...")
    neigh_mat = find_neighbors(x, cut_dist, max_neighs, hard_max = max_neighs_hard_set)
    rows = size(neigh_mat)[1]
    fill!(integrator.p[2][end], 0) # = zeros(size(integrator.p[2][end]))
    integrator.p[2][end][1:rows,:] = neigh_mat
    print("\tDone\n")
end
