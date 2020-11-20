# exports
export minimize

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
            check_units(1energy_unit(pot), UNITS.energy)
            check_units(1force_unit(pot), UNITS.force)
        end
    end

    # >> remove units
    interatomic_potentials = [copy_pot(pot) for pot in interatomic_potentials]
    u0, v0, mass, Δτ = [1ustrip(i) for i in [u0, v0, mass, Δτ]]
    # << remove units

    others = NeighParams(max_neighs, max_neighs_hard_set, reneighboring_every)

    MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids = a_ids, m_ids = m_ids, Δτ = Δτ, save_every = save_every, thermo_save_every = thermo_save_every, others=others)
end

function problem(sim::MDSim, n::Int64, ensemble::Array{<:Ensemble, 1}; callbacks::Array{Function,1}=Function[], verbose::Bool=false)
    tspan = (0.0*sim.Δτ, n*sim.Δτ)
    params = exe_at_start(sim, n, verbose)
    print_thermo_at_start(params, sim, verbose)

    function SOODE(dv, v, u, p, t)
        fill!(dv, 0.0)
        set_acceleration!(dv, v, u, params, sim)
        for ens in ensemble
            ddu!(dv, v, u, params, t, ens)
        end
    end
    cbs = [i(params, ensemble) for i in callbacks]
    p = Float64[]
    prob = MDBase.SecondOrderODEProblem(SOODE, sim.v0, sim.u0, tspan, p, callback=    mdcallbackset(params, ensemble, cbs))
    return prob, sim.Δτ, tspan[1]:sim.save_every*sim.Δτ:tspan[2], params
end

function simulate(n::Int64, sim::MDSim, ensemble::Array{<:Ensemble, 1}; callbacks=Function[], verbose::Bool=false, kwargs...)
    prob, dt, saveat, params = MDSimulator.problem(sim, n, ensemble, callbacks=callbacks, verbose=verbose)
    sol = MDBase.solve(prob, MDBase.VelocityVerlet(), dt=sim.Δτ, saveat=saveat, kwargs...)
    params.M.step = 0
    return sol, params
end

function minimize(sim::MDSim; kwargs...)
    println("Minimizing...")
    parameters = exe_at_start(sim, 1, true)
    loss = (x) -> get_potential_energy(zeros(size(x)...), x, parameters, sim)
    res = optimize(loss, sim.u0, LBFGS(), kwargs...)
    sim.u0 .= res.minimizer
    return sim
end


function MDBase.integrator(sim::MDSim, n::Int64, ensemble::Array{<:Ensemble, 1}; cbs=Function[], verbose::Bool=false)
    prob, dt, saveat, params = problem(sim, n, ensemble, callbacks=cbs, verbose=verbose)
    return MDBase.init(prob, MDBase.VelocityVerlet(), dt=dt), params
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

    cut_dist = ustrip(maximum([cutoff(pot) for pot in sim.interatomic_potentials])*1.1)

    if !(isapprox(cut_dist, 0.0))
        if verbose
            println("Cutoff for reneighboring is set to be $(round(cut_dist*100)/100) $(UNITS.distance).")
        end
        neighs = find_neighbors(sim.u0, cut_dist, sim.boundary_condition, sim.others.max_neighs, hard_max=sim.others.max_neighs_hard_set)
    else
        neighs = zeros(Int64, (0,0))
    end

    nn = fld(n, sim.thermo_save_every) + 1
    thermo_vals = [zeros(nn), zeros(nn), zeros(nn)]

    others = Neighbours(neighs, cut_dist, thermo_vals)

    N = Types.I(size(sim.u0, 2))
    mm = fld(n, sim.save_every) + 1
    acc = zeros(Types.F, (size(sim.u0)..., mm))
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
    params.S.others.thermo_vals[3][1] = get_potential_energy(sim.v0, sim.u0, params, sim)
    if verbose
        pe = params.S.others.thermo_vals[3][1]
        ke = params.S.others.thermo_vals[1][1]
        T = params.S.others.thermo_vals[2][1]
        te = ke+pe
        sig = 4
        println("Time: 1×Δτ\tKE: ", round(ke, sigdigits=sig), "\tTemp: ", round(T, sigdigits=sig), "\tPE: ", round(pe, sigdigits=sig), "\tTE: ", round(te, sigdigits=sig))
    end
end
