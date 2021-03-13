# exports
export potential_force, potential_energy, energy_unit, force_unit

######################################################
# Lennard Jones Potential                            #
######################################################

export LennardJonesParameters

struct LennardJonesParameters{eType, lType, expType <: Real} <: PairPotential
    ϵ::eType
    σ::lType
    R::lType
    exp1::expType
    exp2::expType
    equillibrium::lType
    energy_shift::eType
    ids::Array{Int64,1}
end

function LennardJonesParameters(ϵ::Array{T1,1}, σ::Array{T2,1}, R::Array{T3,1}, exp1::T, exp2::T; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number, T}

	if (Unitful.dimension(first(σ)) == Unitful.dimension(first(R))) && (Unitful.dimension(first(R)) != Unitful.NoDims)
        σ, R = uconvert.(unit(first(σ)), promote(σ, R))
    else
        nothing
    end

	σ = MDSimulator.SymMat(σ)
	ϵ = MDSimulator.SymMat(ϵ)
	R = MDSimulator.SymMat(R)

    n = R.s
    if ids==nothing
        ids = collect(1:n)
    end
    m = maximum(ids)
    σ_ = SymMat([0.0 for i in 1:m for j in 1:i]*1unit(σ))
    ϵ_ = SymMat([0.0 for i in 1:m for j in 1:i]*1unit(ϵ))
    R_ = SymMat([0.0 for i in 1:m for j in 1:i]*1unit(R))
    for ii in 1:n
        for jj in 1:n
            i = ids[ii]
            j = ids[jj]
            σ_[i,j] = σ[ii,jj]
            ϵ_[i,j] = ϵ[ii,jj]
            R_[i,j] = R[ii,jj]
        end
    end

    A = (σ_/R_)^exp1
    B = (σ_/R_)^exp2
    equi_dist = σ_/((exp2/exp1)^(1/(exp1-exp2)))

    LennardJonesParameters(ϵ_, σ_, R_, exp1, exp2, equi_dist, (B-A)*4ϵ_, ids)
end

function LennardJonesParameters(ϵ::Array{T1,1}, σ::Array{T2,1}, R::Array{T3,1}; ids=nothing) where {T1<:Number, T2<:Number, T3<:Number}
    LennardJonesParameters(ϵ, σ, R, 12.0, 6.0, ids=ids)
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number, exp1::T, exp2::T; ids=nothing) where T
    LennardJonesParameters([ϵ], [σ], [R], exp1, exp2, ids=ids)
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number; ids=nothing)
    LennardJonesParameters([ϵ], [σ], [R], ids=ids)
end

function LennardJonesParameters(;ϵ=1.0u"kJ/mol", σ= 100.0u"pm", R=250.0u"pm", ids=nothing)
    LennardJonesParameters(ϵ, σ, R, ids=ids)
end

function Base.similar(lj::LennardJonesParameters, args)
    LennardJonesParameters(args...)
end

@inline function potential_energy(r::F, pot::LennardJonesParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        σ_r = pot.σ[i,j]/r
        A = (σ_r)^pot.exp1
        B = (σ_r)^pot.exp2
        A_ = (pot.σ[i,j]/pot.R[i,j])^pot.exp1
        B_ = (pot.σ[i,j]/pot.R[i,j])^pot.exp2
        return (A-B)*(4pot.ϵ[i,j])
    else
        return 0.0
    end
end

function potential_force(r::F, pot::LennardJonesParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if i<=pot.R.s && j<=pot.R.s && r<pot.R[i,j]
        σ_r = pot.σ[i,j]/r
        A = (σ_r)^pot.exp1
        B = (σ_r)^pot.exp2
        return  (pot.exp1*A-pot.exp2*B)*(4pot.ϵ[i,j])/r
    else
        return 0.0
    end
end

energy_unit(pot::LennardJonesParameters) = unit(pot.ϵ[1])
force_unit(pot::LennardJonesParameters) = unit(pot.ϵ[1]/pot.R[1])
potential_expression(pot::LennardJonesParameters) = "V(r) = 4ϵ[(σ/r)ᵅ - (σ/r)ᵝ]"

######################################################
# Lennard Jones Potential                            #
######################################################
