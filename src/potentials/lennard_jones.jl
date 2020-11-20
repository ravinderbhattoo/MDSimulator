# exports
export potential_force, potential_energy, energy_unit, force_unit

######################################################
# Lennard Jones Potential                            #
######################################################

export LennardJonesParameters

struct LennardJonesParameters{eType <: Number, lType <: Number, expType <: Real} <: PairPotential
    ϵ::Array{eType,2}
    σ::Array{lType,2}
    R::Array{lType,2}
    exp1::expType
    exp2::expType
    equillibrium::Array{lType,2}
    energy_shift::Array{eType,2}
end

function LennardJonesParameters(ϵ::Array{T1,2}, σ::Array{T2,2}, R::Array{T3,2}, exp1::T, exp2::T) where {T1<:Number, T2<:Number, T3<:Number, T}
    A = (σ./R).^exp1
    B = (σ./R).^exp2
    equi_dist = σ./(exp2./exp1).^(1/(exp1-exp2))
    σ, R, equi_dist = uconvert.(unit(σ[1]), promote(σ, R, equi_dist))
    LennardJonesParameters(ϵ, σ, R, exp1, exp2, equi_dist, (B.-A).*4ϵ)
end

function LennardJonesParameters(ϵ::Array{T1,1}, σ::Array{T2,1}, R::Array{T3,1}, exp1::T, exp2::T) where {T1<:Number, T2<:Number, T3<:Number, T}
    ONE = ones(length(σ))
    ϵ_ = sqrt.(ϵ*ϵ')
    σ_ = (σ*ONE' + ONE*σ' ) / 2
    R_ = (R*ONE' + ONE*R' ) / 2
    LennardJonesParameters(ϵ_, σ_, R_, exp1, exp2)
end

function LennardJonesParameters(ϵ::Array{T1,1}, σ::Array{T2,1}, R::Array{T3,1}) where {T1<:Number, T2<:Number, T3<:Number}
    LennardJonesParameters(ϵ, σ, R, 12.0, 6.0)
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number, exp1::T, exp2::T) where T
    LennardJonesParameters([ϵ], [σ], [R], exp1, exp2)
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number)
    LennardJonesParameters([ϵ], [σ], [R])
end

function LennardJonesParameters()
    LennardJonesParameters(1.0u"kJ/mol", 1.0u"pm", 2.5u"pm")
end

function Base.similar(lj::LennardJonesParameters, args)
    LennardJonesParameters(args...)
end

@inline function potential_energy(r::F, pot::LennardJonesParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if r<pot.R[i,j]
        σ_r = pot.σ[i,j]/r
        return ((σ_r)^pot.exp1-(σ_r)^pot.exp2)*(4pot.ϵ[i,j]) + pot.energy_shift[i,j]
    else
        return 0.0
    end
end

function potential_force(r::F, pot::LennardJonesParameters, pair::Tuple{Int64,Int64}) where F <: Number
    i,j = pair
    if r<pot.R[i,j]
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

function Base.show(stream::IO, pp::LennardJonesParameters)
    println(stream, "Lennard-Jones:")
    println(stream, "\tV(r) = 4ϵ[(σ/r)ᵅ - (σ/r)ᵝ] where" )
    function preetify(symb)
        obj = getproperty(pp, symb)
        print(stream, "\t",symb,"(",unit(obj),"):\t")
        show(stream, 1ustrip(obj))
        println(stream)
    end
    for i in fieldnames(typeof(pp))
        preetify(i)
    end
    print(stream, "\tEnergy shift:\t"); show(stream, pp.energy_shift./pp.ϵ); print(stream, "×ϵ"); println(stream)
end

function Base.:+(a::LennardJonesParameters, b::LennardJonesParameters)
    ϵ = [diag(a.ϵ)...; diag(b.ϵ)...]
    σ = [diag(a.σ)...; diag(b.σ)...]
    R = [diag(a.R)...; diag(b.R)...]
    exp1 = a.exp1
    exp2 = a.exp2
    LennardJonesParameters(ϵ, σ, R, exp1, exp2)
end

######################################################
# Lennard Jones Potential                            #
######################################################
