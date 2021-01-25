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
            thermo_vals = params.S.others.thermo_vals
            n = fld(params.M.step, params.S.sim.thermo_save_every) + 1
            thermo_vals.pe[n] = get_potential_energy(v, x, params)
            thermo_vals.ke[n] = params.M.ke
            thermo_vals.temp[n] = params.M.Temperature
            for (key, val) in thermo_vals.custom
                val[n] = params.S.cthermo[key](v, x, params)
            end
            if params.S.verbose
                thermo_print(thermo_vals, n)
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
