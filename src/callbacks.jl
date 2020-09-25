export cb_reneighboring, cb_savethermo

############################################################
# >>>>> Save thermo                                        #
############################################################
function savethermo_condition_f(u,t,integrator)
    # optimized for speed
    ensemble, params = integrator.p.ensemble, integrator.p.params
    (params.M.step)%params.S.sim.thermo_save_every==0
end
function savethermo_affect_f!(integrator)
    # optimized for speed
    ensemble, params = integrator.p.ensemble, integrator.p.params
    x = integrator.u.x[2]
    N = params.S.N
    n = fld(params.M.step, params.S.sim.thermo_save_every)
    params.S.others.thermo_vals[3][n] = get_potential_energy(x, params, params.S.sim)
    params.S.others.thermo_vals[1][n] = params.M.ke
    params.S.others.thermo_vals[2][n] = params.M.Temperature
    if params.S.verbose
        println("Time: ", params.M.step, "×Δτ\tKE: ", params.S.others.thermo_vals[1][n], "\tTemp: ", params.S.others.thermo_vals[2][n], "\tPE: ", params.S.others.thermo_vals[3][n])
    end
    nothing
end
cb_savethermo = MDBase.DiscreteCallback(savethermo_condition_f, savethermo_affect_f!, save_positions=(false,false))

############################################################
# <<<<< Save thermo                                        #
############################################################

############################################################
# >>>>> Reneighboring                                      #
############################################################
function reneighboring_condition_f(u,t,integrator)
    # optimized for speed
    ensemble, params = integrator.p.ensemble, integrator.p.params
    params.M.step % params.S.sim.others.reneighboring_every==0
end

function reneighboring_affect_f!(integrator)
    # optimized for speed
    ensemble, params = integrator.p.ensemble, integrator.p.params
    N = params.S.N
    x = integrator.u.x[2]
    reneighboring!(x, integrator, params)
end

cb_reneighboring = MDBase.DiscreteCallback(reneighboring_condition_f, reneighboring_affect_f!, save_positions=(false,false))
############################################################
# <<<<< Reneighboring                                      #
############################################################
