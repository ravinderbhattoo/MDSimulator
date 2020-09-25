# exports

function Base.show(stream::IO, sim::T) where T <: MDSim
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
    println(stream, "Max. neighbors: $(max(sim.others.max_neighs,sim.others.max_neighs_hard_set))")
    println(stream, "Reneighboring every $(sim.others.reneighboring_every) steps")
    println(stream, "Output:")
    println(stream, "\tLocal (Per atom) save every:\t$(sim.save_every)×Δτ")
    println(stream, "\tthermo (system average) save every:\t$(sim.thermo_save_every)×Δτ")
    println(stream, "Interactions:")
    for (ind, pot) in enumerate(sim.interatomic_potentials)
        print(stream, "$ind. $(pot)")
    end
    println(stream, "<<<<<<<<<<< Simulation Parameters")
end

struct NeighParams{I <: IntType} <: AbstractMDParams
    max_neighs::I
    max_neighs_hard_set::I
    reneighboring_every::I
end

function MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, verify_units::Bool; a_ids = nothing, m_ids = nothing, Δτ = 1.0UNITS.time, save_every = 1000, thermo_save_every = 100, max_neighs = 100, max_neighs_hard_set = 100, reneighboring_every = 100)

    if verify_units
        check_units(u0, UNITS.distance)
        check_units(v0, UNITS.velocity)
        check_units(mass, UNITS.mass)
        check_units(Δτ, UNITS.time)

        for pot in interatomic_potentials
            check_units(potential_energy(0.1UNITS.distance, pot), UNITS.energy)
            check_units(potential_force(0.1UNITS.distance, pot), UNITS.force)
        end
    end

    # >> remove units
    interatomic_potentials_ = [copy_pot(pot) for pot in interatomic_potentials]
    u0, v0, mass, Δτ = [1ustrip(i) for i in [u0, v0, mass, Δτ]]

    # << remove units
    others = NeighParams(max_neighs, max_neighs_hard_set, reneighboring_every)

    print(typeof(interatomic_potentials_[1]))

    MDSim(u0, v0, mass, interatomic_potentials_, boundary_condition; a_ids = a_ids, m_ids = m_ids, Δτ = Δτ, save_every = save_every, thermo_save_every = thermo_save_every, others=others)
end

function simulate(sim::MDSim, n::Int64, ensemble::Array{<:Ensemble, 1}; verbose::Bool=false)
    tspan = (0.0*sim.Δτ, n*sim.Δτ)
    params = exe_at_start(sim, n, verbose)
    print_thermo_at_start(params, sim, verbose)
    p = MDBase.IntgParams(ensemble, params)
    prob = MDBase.SecondOrderODEProblem(MDBase.SOODE, sim.v0, sim.u0, tspan, p, callback=mdcallbackset((cb_savethermo, cb_reneighboring)))
    MDBase.solve(prob, MDBase.VelocityVerlet(), dt=sim.Δτ, saveat=tspan[1]:sim.save_every*sim.Δτ:tspan[2])
end



function MDBase.integrator(sim::MDSim, n::Int64, ensemble::Array{<:Ensemble, 1}; verbose::Bool=false)
    tspan = (0.0*sim.Δτ, n*sim.Δτ)
    params = exe_at_start(sim, n, verbose)
    print_thermo_at_start(params, sim, verbose)
    p = MDBase.IntgParams(ensemble, params)
    prob = MDBase.SecondOrderODEProblem(MDBase.SOODE, sim.v0, sim.u0, tspan, p, callback=mdcallbackset((cb_savethermo, cb_reneighboring)))
    MDBase.init(prob, MDBase.VelocityVerlet(), dt=sim.Δτ)
end


struct Neighbours{F <: FloatType, I <: IntType} <: AbstractMDParams
    neighs::Array{I, 2}
    cut_dist::F
    thermo_vals::Array{Array{F,1} ,1}
end

struct StaticParams{SimObj <: AbstractSimObj, I <: IntType, F <: FloatType} <: AbstractMDParams
    @SParams_generic_fields
    others::Neighbours
    verbose::Bool
end


function exe_at_start(sim::MDSim, n::Types.I, verbose::Bool)
    ux = @view sim.u0[1, :]
    uy = @view sim.u0[2, :]
    uz = @view sim.u0[3, :]
    vx = @view sim.v0[1, :]
    vy = @view sim.v0[2, :]
    vz = @view sim.v0[3, :]
    apply_simulation_bc!(ux, uy, uz, vx, vy, vz, sim.boundary_condition)

    cut_dist = maximum([cutoff(pot) for pot in sim.interatomic_potentials])*1.1

    if !(isapprox(cut_dist, 0.0))
        if verbose
            println("Cutoff for reneighboring is set to be $(round(cut_dist*100)/100) $(UNITS.distance).")
        end
        neighs = find_neighbors(sim.u0, cut_dist, sim.boundary_condition, sim.others.max_neighs, hard_max=sim.others.max_neighs_hard_set)
    else
        neighs = nothing
    end

    nn = fld(n, sim.thermo_save_every) + 1
    thermo_vals = [zeros(nn), zeros(nn), zeros(nn)]

    others = Neighbours(neighs, cut_dist, thermo_vals)

    N = Types.I(size(sim.u0, 2))
    acc = zeros(Types.F, (size(sim.u0)..., n+1))
    mf_acc = unit_factor(UNITS.acceleration, UNITS.force/UNITS.mass)
    kb = 1ustrip(CONSTANTS.kb)


    S = StaticParams(sim, N, kb, mf_acc, acc, others, verbose)
    M = MParams(Types.I(0), Types.F(0), Types.F(0),)

    M.ke = 0.5sum( @. sim.mass*(vx^2 + vy^2 + vz^2) )
    S.others.thermo_vals[1][1] = M.ke
    M.Temperature = 2S.others.thermo_vals[1][1]/(3N*kb)
    S.others.thermo_vals[2][1] = M.Temperature
    return MDParams(S, M)
end

function print_thermo_at_start(params, sim::MDSim, verbose::Bool)
    params.S.others.thermo_vals[3][1] = get_potential_energy(sim.u0, params, sim)
    if verbose
        println("Time: 1×Δτ\tKE: ", params.S.others.thermo_vals[1][1], "\tTemp: ", params.S.others.thermo_vals[2][1], "\tPE: ", params.S.others.thermo_vals[3][1])
    end
end
