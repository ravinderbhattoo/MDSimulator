# exports
export potential_force, potential_energy, energy_unit, force_unit

######################################################
# Electrostatic Potential                            #
######################################################

export ElectrostaticPotentialParameters

# No doc require
struct ElectrostaticPotentialParameters{aType<:SymMat, bType, cType<:SymMat} <: PairPotential
    q₁q₂::aType
    ϵ0::bType
    R::cType
    ids::Array{Int64,1}
end
	
"""
    ElectrostaticPotentialParameters(q::Array{T1,1}, R::Array{T3,1}, ϵ0::T2; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number}

# Arguments
- `q`: Array of charge values i.e. [q₁, q₂, ...]
- `R`: Cutoff values ([R₁₁, R₂₂, ...]) 
- `ϵ0`: ϵ0 

# Keyword Arguments
- `ids`: Ids of the atoms (default: nothing) ([1,2])

# Output
- `struct object`

# Example(s)

ElectrostaticPotentialParameters([1,1], [1,1], [1,1]; ids=[1,2])

"""
function ElectrostaticPotentialParameters(q::Array{T1,1}, R::Array{T3,1}, ϵ0::T2; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number}
    n = length(q)
    if ids==nothing
        ids = collect(1:n)
    end
    m = maximum(ids)
    q_ = SymMat([0.0 for i in 1:m for j in 1:i])
    R_ = SymMat([0.0 for i in 1:m for j in 1:i])
    for ii in 1:n
        for jj in 1:n
            i = ids[ii]
            j = ids[jj]
            q_[i,j] = q[ii]*q[jj] 
            R_[i,j] = (R[ii]+R[jj])/2
        end
    end
        
    ElectrostaticPotentialParameters(q_, ϵ0, R_, ids)
end

function ElectrostaticPotentialParameters(q::Number, R::Number, ϵ0::Number; ids=nothing)
    ElectrostaticPotentialParameters([q], [R], ϵ0, ids=ids)
end

function ElectrostaticPotentialParameters()
    ElectrostaticPotentialParameters(1.0u"C", 1.0u"pm", 1.0u"C*C*mol/kJ/pm", ids=nothing)
end

function Base.similar(lj::ElectrostaticPotentialParameters, args)
    ElectrostaticPotentialParameters(args...)
end

@inline function potential_energy(r::F, pot::ElectrostaticPotentialParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        return pot.q₁q₂[i,j]/(4pi*pot.ϵ0*r)
    else
        return 0.0
    end
end

function potential_force(r::F, pot::ElectrostaticPotentialParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        return  -pot.q₁q₂[i,j]/(4pi*pot.ϵ0*r^2)
    else
        return 0.0
    end
end


energy_unit(pot::ElectrostaticPotentialParameters) = unit(pot.q₁q₂[1,1]/pot.ϵ0/pot.R[1,1])
force_unit(pot::ElectrostaticPotentialParameters) = unit(1*energy_unit(pot)/pot.R[1])
potential_expression(pot::ElectrostaticPotentialParameters) = "V(r) = q₁q₂/4πϵ₀r" 

######################################################
# Electrostatic Potential                            #
######################################################
