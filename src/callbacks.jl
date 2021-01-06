export cb_reneighboring_f, cb_savethermo_f

############################################################
# >>>>> Save thermo                                        #
############################################################
function savethermo_condition_f(params, ensemble)
    return (u,t,integrator) -> begin
            # optimized for speed
            (params.M.step)%params.S.sim.thermo_save_every==0
        end
end
function savethermo_affect_f!(params, ensemble)
    return (integrator) -> begin
            # optimized for speed
            x = integrator.u.x[2]
            v = integrator.u.x[1]
            N = params.S.N
            n = fld(params.M.step, params.S.sim.thermo_save_every) + 1
            params.S.others.thermo_vals[3][n] = get_potential_energy(v, x, params)
            params.S.others.thermo_vals[1][n] = params.M.ke
            params.S.others.thermo_vals[2][n] = params.M.Temperature
            if params.S.verbose
                println("Time: ", params.M.step, "×Δτ\tKE: ", round(params.S.others.thermo_vals[1][n], sigdigits=3), "\tTemp: ", round(params.S.others.thermo_vals[2][n], sigdigits=3), "\tPE: ", round(params.S.others.thermo_vals[3][n], sigdigits=3), "\tTE: ", round(params.S.others.thermo_vals[3][n]+params.S.others.thermo_vals[1][n], sigdigits=3))
            end
            nothing
        end
end
cb_savethermo_f(params, ensemble) = MDBase.DiscreteCallback(savethermo_condition_f(params, ensemble), savethermo_affect_f!(params, ensemble), save_positions=(false,false))

############################################################
# <<<<< Save thermo                                        #
############################################################

############################################################
# >>>>> Reneighboring                                      #
############################################################
function reneighboring_condition_f(params, ensemble)
    return (u,t,integrator) -> begin
            # optimized for speed
            params.M.step % params.S.sim.others.reneighboring_every==0
        end
end

function reneighboring_affect_f!(params, ensemble)
    return (integrator) -> begin
            # optimized for speed
            x = integrator.u.x[2]
            reneighboring!(x, params)
        end
end

cb_reneighboring_f(params, ensemble) = MDBase.DiscreteCallback(reneighboring_condition_f(params, ensemble), reneighboring_affect_f!(params, ensemble), save_positions=(false,false))
############################################################
# <<<<< Reneighboring                                      #
############################################################
