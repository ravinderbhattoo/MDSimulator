export NVE, NVT

"""
  function NVE()

Create an NVE ensemble.
"""
function NVE()
    ENS([],[])
end

#######################################################
# >>>>>>> NVT                                         #
#######################################################

"""
  function NVT(;temperature=300.0)

Create an NVT ensemble with default temperature 300 K with Andersen thermostat.
"""
function NVT(;temperature=300.0)
    ENS([AndersenThermostat(temperature)])
end

"""
  function NVT(thermostats::Array{<:Union{CThermostat, EThermostat},1})

Create an NVT ensemble with given array of thermostats.
"""
function NVT(thermostats::Array{<:Union{CThermostat, EThermostat},1})
    ENS(thermostats)
end
