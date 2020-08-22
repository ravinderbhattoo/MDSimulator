# exports
export SimulationEnvironment
export ddu, ddu!
export NVE

abstract type SimulationEnvironment end


#######################################################
# >>>>>>> NVE                                         #
#######################################################

struct NVE <: SimulationEnvironment
end

function ddu(env::NVE)
    0.0
end

function ddu!(ddu, env::NVE)
    ddu
end

#######################################################
# <<<<<<< NVE                                         #
#######################################################

#######################################################
# >>>>>>> NVT                                         #
#######################################################

struct NVT <: SimulationEnvironment
    thermostat :: Thermostat
    params :: Vector{Float64}
end

function ddu(env::NVT)
    ddu(NVT.thermostat, NVT.params...)
end

function ddu!(ddu, env::NVT)
    @. ddu + ddu(NVT.thermostat, NVT.params...)
end

#######################################################
# <<<<<<< NVT                                         #
#######################################################
