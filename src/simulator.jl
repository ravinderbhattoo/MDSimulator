# exports
export minimize

function Base.show(stream::IO, sim::T) where T <: MDSim
    println(stream, "Simulation Parameters:")
    println(stream, "Units:")
    for f in fieldnames(MDSimulator._UNITS)
        print(stream, "\t$(f):")
        println(stream, "\t$(getproperty(UNITS,f))")
    end
    println(stream, "\nBoundary condition:")
    println(stream, "\t x-direction :\t$(sim.boundary_condition.X)")
    println(stream, "\t y-direction :\t$(sim.boundary_condition.Y)")
    println(stream, "\t z-direction :\t$(sim.boundary_condition.Z)")
    println(stream, "Output:")
        println(stream, "\tLocal (Per atom) save every:\t$(sim.save_every)×Δτ")
        println(stream, "\tthermo (system average) save every:\t$(sim.thermo_save_every)×Δτ")
    println(stream, "\nInteractions:")
    for (ind, pots) in enumerate(sim.interatomic_potentials)
        println(stream, "$ind ", supertype(typeof(pots[1])), "\t→\t")
        for (ind2, pot) in enumerate(pots)
            println(stream, "\n", "$ind.$ind2 ", pot)
        end
    end

    println(stream, "\nOthers:")
    println(stream, "\tTime step (Δτ): $(sim.Δτ)")
    println(stream, "\tMax. neighbors: $(max(sim.others.max_neighs,sim.others.max_neighs_hard_set))")
    println(stream, "\tReneighboring every: $(sim.others.reneighboring_every) steps")

    println(stream, "")
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

function problem(n::Int64, sim::MDSim, ensemble::Array{<:Ensemble, 1}; callbacks::Array{Function,1}=Function[], verbose::Bool=false, cthermo=Dict{String, Function}())
    tspan = (0.0*sim.Δτ, n*sim.Δτ)
    params = exe_at_start(n, sim, verbose, cthermo)
    print_thermo_at_start(params, verbose)

    function SOODE(dv, v, u, p, t)
        fill!(dv, 0.0)
        set_acceleration!(dv, v, u, params)
        for ens in ensemble
            ddu!(dv, v, u, params, t, ens)
        end
        dv .*= params.S.mf_acc
    end
    cbs = [i(params, ensemble) for i in callbacks]
    p = Float64[]
    prob = MDBase.SecondOrderODEProblem(SOODE, sim.v0, sim.u0, tspan, p, callback=    mdcallbackset(params, ensemble, cbs))
    return prob, sim.Δτ, tspan[1]:sim.save_every*sim.Δτ:tspan[2], params
end

function simulate(n::Int64, sim::MDSim, ensemble::Array{<:Ensemble, 1}; callbacks=Function[], verbose::Bool=false, cthermo=Dict{String, Function}(), kwargs...)
    prob, dt, saveat, params = MDSimulator.problem(n, sim, ensemble, callbacks=callbacks, verbose=verbose, cthermo=cthermo)
    sol = MDBase.solve(prob, MDBase.VelocityVerlet(), dt=sim.Δτ, saveat=saveat, kwargs...)
    params.M.step = 0
    (sol.retcode)
    return sol, params
end

function minimize(sim::MDSim; kwargs...)
    println("Minimizing...")
    parameters = exe_at_start(1, sim, true)
    loss = (x) -> get_potential_energy(zeros(size(x)...), x, parameters, sim)
    res = optimize(loss, sim.u0, LBFGS(), kwargs...)
    sim.u0 .= res.minimizer
    return sim
end


function MDBase.integrator(n::Int64, sim::MDSim, ensemble::Array{<:Ensemble, 1}; callbacks=Function[], verbose::Bool=false, cthermo=Dict{String, Function}())
    prob, dt, saveat, params = problem(n, sim, ensemble, callbacks=callbacks, verbose=verbose, cthermo=cthermo)
    return MDBase.init(prob, MDBase.VelocityVerlet(), dt=dt), params
end

struct ThermoVals
    custom::Dict{String,Array{Float64,1}}
    temp::Array{Float64, 1}
    pe::Array{Float64, 1}
    ke::Array{Float64, 1}
end

struct Neighbours{F <: FloatType, I <: IntType} <: AbstractMDParams
    neighs::Array{I, 2}
    cut_dist::F
    thermo_vals::ThermoVals
end

struct StaticParams{I <: IntType, SimObj <: AbstractSimObj, F <: FloatType} <: AbstractMDParams
    @SParams_generic_fields
    others::Neighbours
    verbose::Bool
    cthermo::Dict{String, Fun} where Fun<:Function
end


function exe_at_start(n::Types.I, sim::MDSim, verbose::Bool, cthermo)
    ux = @view sim.u0[1, :]
    uy = @view sim.u0[2, :]
    uz = @view sim.u0[3, :]
    vx = @view sim.v0[1, :]
    vy = @view sim.v0[2, :]
    vz = @view sim.v0[3, :]
    apply_simulation_bc!(ux, uy, uz, vx, vy, vz, sim.boundary_condition)

    cut_dist = ustrip(maximum([cutoff(pot) for pots in sim.interatomic_potentials for pot in pots])*1.1)

    if !(isapprox(cut_dist, 0.0))
        if verbose
            println("Cutoff for reneighboring is set to be $(round(cut_dist*100)/100) $(UNITS.distance).")
        end
        neighs = find_neighbors(sim.u0, cut_dist, sim.boundary_condition, sim.others.max_neighs, hard_max=sim.others.max_neighs_hard_set)
    else
        neighs = zeros(Int64, (0,0))
    end

    nn = fld(n, sim.thermo_save_every) + 1
    custom = Dict{String, Array{Float64,1}}()
    for key in keys(cthermo)
        custom[key] = zeros(nn)
    end
    thermo_vals = ThermoVals(custom, zeros(nn), zeros(nn), zeros(nn))

    others = Neighbours(neighs, cut_dist, thermo_vals)

    N = Types.I(size(sim.u0, 2))
    mm = fld(n, sim.save_every) + 1
    acc = zeros(Types.F, (size(sim.u0)..., mm))
    mf_acc = unit_factor(UNITS.acceleration, UNITS.force/UNITS.mass)
    kb = 1ustrip(CONSTANTS.kb)

    S = StaticParams(N, sim, kb, mf_acc, acc, others, verbose, cthermo)
    M = MParams(Types.I(0), Types.F(0), Types.F(0),)

    M.ke = 0.5sum( @. sim.mass*(vx^2 + vy^2 + vz^2) )
    S.others.thermo_vals.ke[1] = M.ke
    M.Temperature = 2S.others.thermo_vals.ke[1]/(3N*kb)
    S.others.thermo_vals.temp[1] = M.Temperature
    return MDParams(S, M)
end

function print_thermo_at_start(params, verbose::Bool)
    thermo_vals = params.S.others.thermo_vals
    thermo_vals.pe[1] = get_potential_energy(params.S.sim.v0, params.S.sim.u0, params)
    for (k, v) in thermo_vals.custom
        v[1] = params.S.cthermo[k](params.S.sim.v0, params.S.sim.u0, params)
    end

    if verbose
        thermo_headers()
        thermo_print(params.S.others.thermo_vals, 1)
    end
end
