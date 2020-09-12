# exports
export MDSim, simulate, sim_init

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

function sim_init(sim::MDSim; kwargs...)
    args = []
    k = keys(kwargs)
    for f in fieldnames(typeof(sim))
        if f in k
            push!(args, kwargs[Symbol(f)])
        else
            push!(args, getproperty(sim,f))
        end
    end
    MDSim(args...)
end

function sim_init(res::MDBase.ODESolution, sim::MDSim)
    data = Array(reshape(Array(res.u[end]),(3,:)))
    N = Int(size(data, 2)/2)
    v0 = data[:, 1:N]
    u0 = data[:, N+1:end]
    sim_init(sim, u0=u0, v0=v0)
end

function MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids = nothing, m_ids = nothing, Δτ = 1.0UNITS.time, save_every = 1000, global_save_every = 100, max_neighs = 100, max_neighs_hard_set = 100, reneighboring_every = 100, check_units=true)

    if check_units
        check_units(u0, UNITS.distance)
        check_units(v0, UNITS.velocity)
        check_units(mass, UNITS.mass)
        check_units(Δτ, UNITS.time)

        for pot in interatomic_potentials
            check_units(potential_energy(0.1UNITS.distance, pot), UNITS.energy)
            check_units(potential_force(0.1UNITS.distance, pot), UNITS.force)
        end
    end

    if a_ids==nothing
        a_ids = ones(Int64, size(u0,2))
        m_ids = copy(a_ids)
    elseif m_ids==nothing
        m_ids = copy(a_ids)
    end

    MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, a_ids, m_ids, Δτ, save_every, global_save_every, max_neighs, max_neighs_hard_set, reneighboring_every)
end


function simulate(Sim::MDSim, n::Int64, ensemble; verbose::Bool=false)
    Δτ = ustrip(Sim.Δτ)
    tspan = (0.0, n*Δτ)

    u0, v0 = 1*ustrip(Sim.u0), 1*ustrip(Sim.v0)
    apply_simulation_bc!(u0, v0, Sim.boundary_condition)
    mass = 1*ustrip(Sim.mass)
    Temp = 1*ustrip(get_temperature(v0, mass))

    cut_dist = ustrip(maximum([cutoff(pot) for pot in Sim.interatomic_potentials]))*1.1

    print("Cutoff for reneighboring is set to be $(round(cut_dist*100)/100) $(UNITS.distance).\n\n")

    neighs = find_neighbors(u0, cut_dist, Sim.boundary_condition, Sim.max_neighs, hard_max=Sim.max_neighs_hard_set)

    nn = Int64(n/Sim.global_save_every) + 1
    global_vals = [WellArray(nn, fill_with=0.0UNITS.energy), WellArray(nn, fill_with=0.0UNITS.temperature), WellArray(nn, fill_with=0.0UNITS.energy)]

    uf = unit_factor(UNITS.acceleration, UNITS.force/UNITS.mass)

    interatomic_potentials = []
    for pot in Sim.interatomic_potentials
        push!(interatomic_potentials, copy_pot(pot))
    end

    dofs = 3*size(u0, 2)

    acceleration = zeros(size(u0)..., length(0.0Δτ:Sim.save_every*Δτ:n*Δτ))

    kb = 1*ustrip(CONSTANTS.kb)

    params = (interatomic_potentials=interatomic_potentials, a_ids=Sim.a_ids, m_ids=Sim.m_ids, mass=mass, boundary_condition=Sim.boundary_condition, neighs=neighs, uf=uf, N=size(u0, 2), cut_dist=cut_dist, max_neighs=Sim.max_neighs, max_neighs_hard_set=Sim.max_neighs_hard_set, reneighboring_every=Sim.reneighboring_every, global_save_every=Sim.global_save_every, global_vals=global_vals, verbose=verbose, dofs=dofs, acc=acceleration, save_every=Sim.save_every, T=@MVector[Temp], nstep=@MVector[0], kb=kb)

    fillit!(params.global_vals[1], unit_convert(UNITS.energy, get_kinetic_energy(v0, mass)UNITS.mass*UNITS.velocity^2))
    fillit!(params.global_vals[2], unit_convert(UNITS.temperature, get_temperature(params.global_vals[1][end], params.dofs)))
    fillit!(params.global_vals[3], get_potential_energy(u0, params)UNITS.energy)

    p = (ensemble, params)

    soode(v, u, p, t) = MD_soode(v, u, p, t)

    prob = SecondOrderODEProblem(soode, v0, u0, tspan, p, callback=mdcallbackset())

    print("Time: 0×Δτ \tKE: $(params.global_vals[1][end]) \tPE: $(params.global_vals[3][end]) \tTemp: $(params.global_vals[2][end]) \n")

    solve(prob, VelocityVerlet(), dt=Δτ, saveat=0.0Δτ:Sim.save_every*Δτ:n*Δτ)
end

function MD_soode(v, u, p, t)
    ensemble, params = p
    dv = get_acceleration(v, u, params)
    for ens in ensemble
        ddu!(dv, v, u, params, t, ens)
    end
    return dv
end
