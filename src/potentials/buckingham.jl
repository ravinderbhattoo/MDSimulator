# exports
export potential_force, potential_energy, energy_unit, force_unit

######################################################
# Buckingham Potential                               #
######################################################

export BuckinghamParameters

# No doc require
struct BuckinghamParameters{AType, BType, CType, lType} <: PairPotential
    A::AType
    B::BType
    C::CType
    R::lType
    equillibrium::lType
    well_energy::AType
    energy_shift::AType
    ids::Array{Int64,1}
end

"""
    BuckinghamParameters(A::T1, B::T2, C::T3, R::T4; ids=nothing) where {T1, T2, T3, T4}

# Arguments
- `A`: Array of pair paramters ([A₁₁, A₁₂, A₂₂])
- `B`: Array of pair paramters ([B₁₁, B₁₂, B₂₂]) 
- `C`: Array of pair paramters ([C₁₁, C₁₂, C₂₂]) 
- `R`: Cutoff values for each pair ([R₁₁, R₁₂, R₂₂])

# Keyword Arguments
- `ids`: Ids of the atoms (default: nothing) ([1,2])

# Output
- `struct object`

# Example(s)

BuckinghamParameters([1,1,1], [1,1,1], [1,1,1], [1,1,1]; ids=[1,2])

"""
function BuckinghamParameters(A::T1, B::T2, C::T3, R::T4; ids=nothing) where {T1, T2, T3, T4}

    n = A.s
    if ids==nothing
        ids = collect(1:n)
    end
    m = maximum(ids)
    A_ = 1unit(A)*SymMat([0.0 for i in 1:m for j in 1:i])
    B_ = 1unit(B)*SymMat([0.0 for i in 1:m for j in 1:i])
    C_ = 1unit(C)*SymMat([0.0 for i in 1:m for j in 1:i])
    R_ = 1unit(R)*SymMat([0.0 for i in 1:m for j in 1:i])
    for ii in 1:n
        for jj in 1:n
            i = ids[ii]
            j = ids[jj]
            A_[i,j] = A[ii,jj] 
            B_[i,j] = B[ii,jj] 
            C_[i,j] = C[ii,jj] 
            R_[i,j] = R[ii,jj]
        end
    end

    energy_shift = A_*exp(-B_*R_) - C_/R_^6
    r = 100R_
    for i in 1:10000
        ∂f_∂r = -A_*B_*exp(-B_*r) + 6C_/r^7
        r -= ∂f_∂r/abs(∂f_∂r)*R_/100
    end
    well_energy = A_*exp(-B_*r) - C_/r^6
    BuckinghamParameters(A_, B_, C_, R_, r, well_energy, -energy_shift, ids)
end

function BuckinghamParameters(A::Array{T1,1}, B::Array{T2,1}, C::Array{T3,1}, R::Array{T4,1}; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}        
    BuckinghamParameters(SymMat(A), SymMat(B), SymMat(C), SymMat(R), ids=ids)
end

function BuckinghamParameters(A::T1, B::T2, C::T3, R::T4; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
    BuckinghamParameters([A], [B], [C], [R], ids=ids)
end

function BuckinghamParameters()
    A = 1.69*1.0e-18u"kJ"*CONSTANTS.Avogadro
    B = 1/27.3u"pm"
    C = 102*1.0e-28u"kJ"*CONSTANTS.Avogadro*(1u"pm")^6
    R = 500u"pm"
    BuckinghamParameters(A, B, C, R, ids=nothing)
end

function potential_energy(r::Number, pot::BuckinghamParameters, pair::Tuple{Int64,Int64})
    i,j  = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        return pot.A[i,j]*(exp(-pot.B[i,j]*r) -  exp(-pot.B[i,j]*pot.R[i,j])) - pot.C[i,j]*(1/r^6 - 1/pot.R[i,j]^6)
    else
        return 0.0
    end
end

function potential_force(r::Number, pot::BuckinghamParameters, pair::Tuple{Int64,Int64})
    i,j = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        return pot.A[i,j]*pot.B[i,j]*exp(-pot.B[i,j]*r)-6pot.C[i,j]/r^7
    else
        return 0.0
    end
end

energy_unit(pot::BuckinghamParameters) = unit(pot.well_energy[1])
force_unit(pot::BuckinghamParameters) = unit(pot.C[1]/pot.R[1]^7)
potential_expression(pot::BuckinghamParameters) = "V(r) = A e⁻ᴮʳ - C r⁻⁶"

######################################################
# Buckingham Potential                               #
######################################################


#
