export mdcallbackset

function mdcallbackset(;list=false)
    cb1 = DiscreteCallback(updates_condition_f, updates_affect_f!, save_positions=(false,false))
    cb2 = DiscreteCallback(saveacc_condition_f, saveacc_affect_f!, save_positions=(false,false))
    cb3 = DiscreteCallback(saveglobal_condition_f, saveglobal_affect_f!, save_positions=(false,false))
    cb4 = DiscreteCallback(reneighboring_condition_f, reneighboring_affect_f!, save_positions=(false,false))
    if !list
        CallbackSet(cb1, cb2, cb3, cb4)
    else
        return cb1, cb2, cb3, cb4
    end
end

function mdcallbackset(cbs::Array{Any,1})
    lisT = mdcallbackset()
    for i in cbs
        push!(lisT, i)
    end
    CallbackSet(lisT...)
end

############################################################
# >>>>> Updates                                            #
############################################################
function updates_condition_f(u,t,integrator)
    true
end
function updates_affect_f!(integrator)
    ensemble, params = integrator.p
    params.nstep[1] += 1
    N = params.N
    data = reshape(Array(integrator.u),(3,:))
    v = data[:,1:N]
    u = data[:,N+1:end]

    # apply boundary conditions
    apply_simulation_bc!(u, v, params.boundary_condition)

    # apply temperature updation required for some tasks
    params.T[1] = 1*ustrip(get_temperature(v, params.mass))

    # apply "apply!" functions from ensembles such as rescaling velocities and positions
    for ens in ensemble
        apply!(v, u, params, integrator.t, ens)
    end

    # updating u and v
    integrator.u[:] = hcat(reshape(v,(1,:)),reshape(u,(1,:)))
    integrator.sol.u[end][:] = integrator.u[:]
end
############################################################
# <<<<< Updates                                            #
############################################################

############################################################
# >>>>> Save acceleration                                  #
############################################################
function saveacc_condition_f(u,t,integrator)
    ensemble, params = integrator.p
    round(integrator.t/integrator.dt)%params.save_every==0
end
function saveacc_affect_f!(integrator)
    ensemble, params = integrator.p
    N = params.N
    n = Int(round(integrator.t/integrator.dt)/params.save_every)
    data = Array(reshape(get_du(integrator),(3,:)))
    params.acc[:,:,n] = data[:,N+1:end]
end
############################################################
# <<<<< Save acceleration                                  #
############################################################

############################################################
# >>>>> Save global                                        #
############################################################
function saveglobal_condition_f(u,t,integrator)
    ensemble, params = integrator.p
    round(integrator.t/integrator.dt)%params.global_save_every==0
end
function saveglobal_affect_f!(integrator)
    ensemble, params = integrator.p
    N = params.N
    data = reshape(Array(integrator.u),(3,:))
    v = data[:,1:N]
    x = data[:,N+1:end]

    fillit!(params.global_vals[1], unit_convert(UNITS.energy, get_kinetic_energy(v, params.mass)UNITS.mass*UNITS.velocity^2))
    fillit!(params.global_vals[2], unit_convert(UNITS.temperature, get_temperature(params.global_vals[1][end], params.dofs)))
    fillit!(params.global_vals[3], get_potential_energy(x, params)UNITS.energy)
    thermo_print(integrator, params)
end
############################################################
# <<<<< Save global                                        #
############################################################

############################################################
# >>>>> Reneighboring                                      #
############################################################
function reneighboring_condition_f(u,t,integrator)
    ensemble, params = integrator.p
    if params.reneighboring_every!=Inf
        isapprox(round(integrator.t/integrator.dt)%params.reneighboring_every,0.0)
    else
        false
    end
end
function reneighboring_affect_f!(integrator)
    params = integrator.p[2]
    N = params.N
    x = reshape(Array(integrator.u),(3,:))[:,N+1:end]
    reneighboring!(x, integrator, params)
end
############################################################
# <<<<< Reneighboring                                      #
############################################################
