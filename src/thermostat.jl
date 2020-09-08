# exports

export CThermostat, EThermostat

export NoseHooverThermostat, AndersenThermostat
export ddu, ddu!, apply!

abstract type CThermostat <: CallBackType end
abstract type EThermostat <: EoMType end

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

@inline function apply!(v, u, params, t, thermo::AndersenThermostat)
    if params.nstep[1]%thermo.ncycle==0 || params.nstep[1]==1
        v .*= sqrt(thermo.T₀/params.T[1])
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

@inline function apply!(v, u, params, t, thermo::MAndersenThermostat)
    if params.nstep[1]%thermo.ncycle==0 || params.nstep[1]==1
        v .*= sqrt((thermo.T₁ + min(thermo.width,thermo.counter)/thermo.width*(thermo.T₂ - thermo.T₁))/params.T[1])
        thermo.counter += thermo.ncycle
    end
end

#######################################################
# <<<<<<< Andersen thermostat                         #
#######################################################

#######################################################
# >>>>>>> Nosé-Hoover thermostat                      #
#######################################################

# http://www.sklogwiki.org/SklogWiki/index.php/Nos%C3%A9-Hoover_thermostat

mutable struct NoseHooverThermostat <: EThermostat
    Q::Real # dimension of energy×(time)²
    X::Int64 # dofs
    ζ::Real # friction
    T0::Real # temperature
end

function NoseHooverThermostat()
    NoseHooverThermostat(1,216*3,1,300)
end

@inline function ddu(v, u, params, t, Obj::NoseHooverThermostat)
    dζ = (sum(params.mass.*sum(v.^2, dims=1)') .- (Obj.X+1)*params.kb*Obj.T0 )/Obj.Q
    Obj.ζ += dζ
    return -Obj.ζ.*v
end

@inline function ddu!(dv, args...)
    dv .+= ddu(args...)
end

#######################################################
# <<<<<<< Nosé-Hoover thermostat                      #
#######################################################
