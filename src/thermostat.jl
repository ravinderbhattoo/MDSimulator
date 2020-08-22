# exports

export Thermostat
export NoseHooverThermostat
export ddu!

abstract type Thermostat end


#######################################################
# >>>>>>> Nosé-Hoover thermostat                      #
#######################################################

# http://www.sklogwiki.org/SklogWiki/index.php/Nos%C3%A9-Hoover_thermostat

struct NoseHooverThermostat <: Thermostat
    Q::Float64 # dimension of energy×(time)²
    X::Int64
    u::Float64 # number of degrees of freedom
end

function ddu!(dv, v, u)
    @. ζ = 1/Q *( mv2 - (X+1)*kB*T )
    @. dv += -ζ*v
end


#######################################################
# <<<<<<< Nosé-Hoover thermostat                      #
#######################################################
