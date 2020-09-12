# exports
export potential_force, potential_energy

######################################################
# Lennard Jones Potential                            #
######################################################

export LennardJonesParameters

struct LennardJonesParameters{eType <: Number, lType <: Number, expType <: Real} <: PairPotential
    ϵ::eType
    σ::lType
    R::lType
    exp1::expType
    exp2::expType
    equillibrium::lType
    energy_shift::eType
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number, exp1::Number, exp2::Number)
    A = (σ/R)^exp1
    B = (σ/R)^exp2
    equi_dist = σ/(exp2/exp1)^(1/(exp1-exp2))
    σ, R, equi_dist = uconvert.(unit(σ), promote(σ, R, equi_dist))
    LennardJonesParameters(ϵ, σ, R, exp1, exp2, equi_dist, (B-A)*4ϵ)
end

function LennardJonesParameters(ϵ::Number, σ::Number, R::Number)
    LennardJonesParameters(ϵ, σ, R, 12.0, 6.0)
end

function LennardJonesParameters()
    LennardJonesParameters(1.0, 1.0, 2.5)
end

function Base.similar(lj::LennardJonesParameters, args)
    LennardJonesParameters(args...)
end

function Base.show(stream::IO, pp::LennardJonesParameters)
    println(stream, "Lennard-Jones:")
    println(stream, "\tV(r) = 4ϵ[(σ/r)ᵅ - (σ/r)ᵝ] where" )
    print(stream, "\tϵ:\t"); show(stream, pp.ϵ); println(stream)
    print(stream, "\tσ:\t"); show(stream, pp.σ); println(stream)
    print(stream, "\tR:\t"); show(stream, pp.R); println(stream)
    print(stream, "\tα:\t"); show(stream, pp.exp1); println(stream)
    print(stream, "\tβ:\t"); show(stream, pp.exp2); println(stream)
    print(stream, "\tEquillibrium position:\t"); show(stream, pp.equillibrium); println(stream)
    print(stream, "\tEnergy shift:\t"); show(stream, pp.energy_shift/pp.ϵ); print(stream, "×ϵ"); println(stream)
end


function potential_energy(r::Number, pot::LennardJonesParameters)
    if unit(r)==unit(pot.σ)
        #
    else
        @warn "Using different units may cause slow code."
    end

    if r<pot.R
        σ_r = pot.σ/r
        return ((σ_r)^pot.exp1-(σ_r)^pot.exp2)*(4pot.ϵ) + pot.energy_shift
    else
        return 0.0
    end
end

function potential_force(r::Number, pot::LennardJonesParameters)
    if unit(r)==unit(pot.σ)
        #
    else
        @warn "Using different units may cause slow code."
    end
    if r<pot.R
        σ_r = pot.σ/r
        A = (σ_r)^pot.exp1
        B = (σ_r)^pot.exp2
        return  (pot.exp1*A-pot.exp2*B)*(4pot.ϵ)/r
    else
        return 0.0
    end
end

######################################################
# Lennard Jones Potential                            #
######################################################
