# exports
export potential_force, potential_energy

######################################################
# NN Potential                                       #
######################################################

export NN

struct NN <: MLPotential
    m::Any
end


function Base.show(stream::IO, pp::NN)
    println(stream, "NN Potential:")
    println(stream, "\tforce = m(x) where" )
    print(stream, "\tm:\t"); show(stream, pp.m); println(stream)
end


function potential_energy(x::Any, pot::NN)
    return 0.0
end

function acceleration(v::Array{<:Number,2}, u::Array{<:Number,2}, pot::NN)
    return pot.m(vcat(u,v))
end

######################################################
# NN Potential                                       #
######################################################
