using DifferentialEquations
using ParticlesTools
using Unitful

# exports
export MDSim, simulate

struct MDSim
    u0
    v0
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

function Base.show(stream::IO, sim::MDSim)
    println(stream, ">>>>>>>>>>> Simulation Parameters")
    println(stream, "Units:")
    for f in fieldnames(MDSimulator._UNITS)
        print(stream, "\t$(f):")
        println(stream, "\t$(getproperty(UNITS,f))")
    end
    println(stream, "Boundary condition:")
    println(stream, "\t x-direction :\t$(sim.boundary_condition.X)")
    println(stream, "\t y-direction :\t$(sim.boundary_condition.Y)")
    println(stream, "\t z-direction :\t$(sim.boundary_condition.Z)")
    println(stream, "Time step (Δτ): $(sim.Δτ)")
    println(stream, "Max. neighbors: $(max(sim.max_neighs,sim.max_neighs_hard_set))")
    println(stream, "Reneighboring every $(sim.reneighboring_every) steps")
    println(stream, "Output:")
    println(stream, "\tLocal (Per atom) save every:\t$(sim.save_every)×Δτ")
    println(stream, "\tGlobal (system average) save every:\t$(sim.global_save_every)×Δτ")
    println(stream, "Interactions:")
    for (ind, pot) in enumerate(sim.interatomic_potentials)
        print(stream, "$ind. $(pot)")
    end
    println(stream, "<<<<<<<<<<< Simulation Parameters")
end



function MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids = nothing, m_ids = nothing, Δτ = 1.0UNITS.time, save_every = 1000, global_save_every = 100, max_neighs = 100, max_neighs_hard_set = 100, reneighboring_every = 100)

    check_units(u0, UNITS.distance)
    check_units(v0, UNITS.velocity)
    check_units(mass, UNITS.mass)
    check_units(Δτ, UNITS.time)

    for pot in interatomic_potentials
        check_units(potential_energy(0.1UNITS.distance, pot), UNITS.energy)
        check_units(potential_force(0.1UNITS.distance, pot), UNITS.force)
    end

    if a_ids==nothing
        a_ids = ones(size(initial_condition)[1][2])
        m_ids = copy(a_ids)
    elseif m_ids==nothing
        m_ids = copy(a_ids)
    end

    MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, a_ids, m_ids, Δτ, save_every, global_save_every, max_neighs, max_neighs_hard_set, reneighboring_every)
end


function MD_soode(v, u, p, t)
    emvironment, params = p
    dv = get_acceleration(u, params)
    # for env in environment
    #     @. dv += ddu(env)
    # end
    return dv
end

function simulate(Sim::MDSim, n::Int64, environment; verbose::Bool=false)
    Δτ = ustrip(Sim.Δτ)
    tspan = (0.0, n*Δτ)

    u0, v0 = 1*ustrip(Sim.u0), 1*ustrip(Sim.v0)
    apply_simulation_bc!(u0, v0, Sim.boundary_condition)
    mass = 1*ustrip(Sim.mass)

    cut_dist = ustrip(maximum([pot.R for pot in Sim.interatomic_potentials]))

    neighs = find_neighbors(u0, cut_dist, Sim.boundary_condition, Sim.max_neighs, hard_max=Sim.max_neighs_hard_set)

    nn = Int64(n/Sim.global_save_every)
    global_vals = [WellArray(nn, fill_with=0.0UNITS.energy), WellArray(nn, fill_with=0.0UNITS.temperature), WellArray(nn, fill_with=0.0UNITS.energy)]

    uf = unit_factor(UNITS.acceleration, UNITS.force/UNITS.mass)

    interatomic_potentials = []
    for pot in Sim.interatomic_potentials
        push!(interatomic_potentials, copy_pot(pot))
    end

    dofs = 3*size(u0, 2)

    params = (interatomic_potentials=interatomic_potentials, a_ids=Sim.a_ids, m_ids=Sim.m_ids, mass=mass, boundary_condition=Sim.boundary_condition, neighs=neighs, uf=uf, N=size(u0, 2), cut_dist=cut_dist, max_neighs=Sim.max_neighs, max_neighs_hard_set=Sim.max_neighs_hard_set, reneighboring_every=Sim.reneighboring_every, global_save_every=Sim.global_save_every, global_vals=global_vals, verbose=verbose, dofs=dofs)

    p = (environment, params)

    soode(v, u, p, t) = MD_soode(v, u, p, t)

    prob = SecondOrderODEProblem(soode, v0, u0, tspan, p, callback=DiscreteCallback(condition_f, affect_f!, save_positions=(false,false)))
    solve(prob, VelocityVerlet(), dt=Δτ, saveat=0.0Δτ:Sim.save_every*Δτ:n*Δτ)
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
        fillit!(params.global_vals[1], unit_convert(UNITS.energy, get_kinetic_energy(v, params.mass)UNITS.mass*UNITS.velocity^2))
        fillit!(params.global_vals[2], unit_convert(UNITS.temperature, get_temperature(params.global_vals[1][end], params.dofs)))
        fillit!(params.global_vals[3], get_potential_energy(x, params)UNITS.energy)
        print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \tKE: $(params.global_vals[1][end]) \tPE: $(params.global_vals[3][end]) \tTemp: $(params.global_vals[2][end]) \n")
    end
    isapprox(round(integrator.t/integrator.dt)%params.reneighboring_every,0.0)
end

function affect_f!(integrator)
    params = integrator.p[2]
    N = params.N
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ \nReneighboring(cutoff = $(params.cut_dist)$(UNITS.distance))...")
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


function copy_pot(pot)
    args = []
    for f in fieldnames(typeof(pot))
        push!(args, ustrip(getproperty(pot,f)))
    end
    similar(pot, args)
end
