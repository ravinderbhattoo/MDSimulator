# exports

export CThermostat, EThermostat

export NoseHooverThermostat, AndersenThermostat
export ddu, ddu!, apply!

abstract type CThermostat <: MDBase.CallBackType end
abstract type EThermostat <: MDBase.EoMType end

#######################################################
# >>>>>>> Andersen thermostat                         #
#######################################################

struct AndersenThermostat <: CThermostat
    T₀::Real # Equilibrium temperature
    ncycle::Int64 # Applied after n steps
end

function AndersenThermostat(temperature::Real; ncycle::Int64=100)
    AndersenThermostat(temperature, ncycle)
end

function Base.show(stream::IO, thermo::AndersenThermostat)
    println(stream, "Anderson thermostat:")
    println(stream, "  Temperature: $(thermo.T₀)")
    println(stream, "  ncycle: $(thermo.ncycle) steps")
end

@inline function MDBase.apply!(v, u, params, t, thermo::AndersenThermostat)
    if params.M.step%thermo.ncycle==0 || params.M.step==1
        v .*= sqrt(thermo.T₀/params.M.Temperature)
    end
end

mutable struct MAndersenThermostat <: CThermostat
    T₁::Real # Equilibrium temperature 1
    T₂::Real # Equilibrium temperature 1
    width::Int64 # will reach t1 to t2 in width steps
    counter::Int64
    ncycle::Int64 # Applied after n steps
end

function AndersenThermostat(t1::Real, t2::Real; width::Int64=10000, ncycle::Int64=100)
    MAndersenThermostat(t1, t2, width, 0, ncycle)
end

function Base.show(stream::IO, thermo::MAndersenThermostat)
    println(stream, "Anderson thermostat (with variable temperature):")
    println(stream, "  Initial Temperature: $(thermo.T₁)")
    println(stream, "  Final Temperature: $(thermo.T₂)")
    println(stream, "  Width: $(thermo.width)")
    println(stream, "  ncycle: $(thermo.ncycle) steps")
end

@inline function MDBase.apply!(v, u, params, t, thermo::MAndersenThermostat)
    if params.M.step%thermo.ncycle==0 || params.M.step==1
        v .*= sqrt((thermo.T₁ + min(thermo.width,thermo.counter)/thermo.width*(thermo.T₂ - thermo.T₁))/params.M.Teperature)
        thermo.counter += thermo.ncycle
    end
end

#######################################################
# <<<<<<< Andersen thermostat                         #
#######################################################
