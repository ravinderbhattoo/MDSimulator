using Unitful

#exports
export UNITS, unit_factor, unit_convert, CONSTANTS

mutable struct _UNITS
    energy
    force
    distance
    time
    velocity
    acceleration
    mass
    temperature
    _UNITS(;energy=u"kJ/mol", distance=u"pm", time=u"fs", mass=u"g/mol", temperature=u"K") = new(energy, energy/distance, distance, time, distance/time, distance/time^2, mass, temperature)
end

UNITS = _UNITS()

function unit_convert(x::Unitful.Units, y::Number)
    a, b = promote(1x, y)
    return uconvert(x, Float64(b/a)*a)
end

function unit_factor(x::Unitful.Units, y::Unitful.Units)
    a, b = promote(1x, 1y)
    return Float64(b/a)
end

Unitful.unit(M::Array{T, N}) where {T, N} = unit(T)

Unitful.ustrip(M::WellArray{T}) where {T} = [ustrip(i) for i in M]
Unitful.ustrip(M::Array{T, N}) where {T, N} = reshape([ustrip(i) for i in M],size(M))

function check_units(a::Any, b::Unitful.Units; name1 = "A", name2 = "B")
    if unit(a) != b
        throw("Unit mismatch. $name1 has unit of $b[$(typeof(b))] while $name2 has unit of $(unit(a))[$(typeof(unit(a)))]")
    end
end


struct _CONSTANTS
    kb
    _CONSTANTS(;kb=8.3144598e-3u"kJ/K/mol") = new(kb)
end

CONSTANTS = _CONSTANTS()

#
