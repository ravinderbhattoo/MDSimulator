abstract type PotentialParameters
end

# exports
export potential_force, potential_energy, PotentialParameters

######################################################
# Lennard Jones Potential                            #
######################################################

export LennardJonesParameters

struct LennardJonesParameters{pType <: Real} <: PotentialParameters
    ϵ::pType
    σ::pType
    R::pType
    exp1::pType
    exp2::pType
    equillibrium::pType
    energy_shift::pType
end

function LennardJonesParameters(ϵ::Real, σ::Real, R::Real, exp1::Real, exp2::Real)
    A = (σ/R)^exp1
    B = (σ/R)^exp2
    LennardJonesParameters(ϵ, σ, R, exp1, exp2, σ/(exp2/exp1)^(1/(exp1-exp2)), (B-A)*4ϵ)
end

function LennardJonesParameters(ϵ::Real, σ::Real, R::Real)
    LennardJonesParameters(ϵ, σ, R, 12.0, 6.0)
end

function LennardJonesParameters()
    LennardJonesParameters(1.0, 1.0, 2.5)
end

function Base.show(stream::IO, pp::LennardJonesParameters)
    println(stream, "Lennard-Jones:")
    print(stream, "\tϵ:"); show(stream, pp.ϵ); println(stream)
    print(stream, "\tσ:"); show(stream, pp.σ); println(stream)
    print(stream, "\tR:"); show(stream, pp.R); println(stream)
    print(stream, "\tExponent_1:"); show(stream, pp.exp1); println(stream)
    print(stream, "\tExponent_2:"); show(stream, pp.exp2); println(stream)
    print(stream, "\tEquillibrium position:"); show(stream, pp.equillibrium); println(stream)
end


function potential_energy(r::Real, pot::LennardJonesParameters)
    if r<pot.R
        σ_r = pot.σ/r
        A = (σ_r)^pot.exp1
        B = (σ_r)^pot.exp2
        return (A-B)*(4pot.ϵ) + pot.energy_shift
    else
        return 0.0
    end
end

function potential_force(r::Real, pot::LennardJonesParameters)
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

######################################################
# Buckingham Potential                               #
######################################################

export BuckinghamParameters

struct BuckinghamParameters{pType <: Real} <: PotentialParameters
    A::pType
    B::pType
    C::pType
    R::pType
    equillibrium::pType
    energy_shift::pType
end

function BuckinghamParameters(A::Real, B::Real, C::Real, R::Real)
    energy_shift = A*exp(-B*R) - C/R^6
    r = R
    for i in 1:100
        r = (6C/(A*B)*exp(B*r))^(1/7)
    end
    BuckinghamParameters(A, B, C, R, r, -energy_shift)
end

function BuckinghamParameters()
    BuckinghamParameters(1.0, 1.0, 0.001, 2.5)
end

function Base.show(stream::IO, pp::BuckinghamParameters)
    println(stream, "Buckingham:")
    print(stream, "\tA:"); show(stream, pp.A); println(stream)
    print(stream, "\tB:"); show(stream, pp.B); println(stream)
    print(stream, "\tC:"); show(stream, pp.C); println(stream)
    print(stream, "\tR:"); show(stream, pp.R); println(stream)
    print(stream, "\tEquillibrium position:"); show(stream, pp.equillibrium); println(stream)
end


function potential_energy(r::Real, pot::BuckinghamParameters)
    if r<pot.R
        return pot.A*exp(-pot.B*r) - pot.C/r^6 + pot.energy_shift
    else
        return 0.0
    end
end

function potential_force(r::Real, pot::BuckinghamParameters)
    if r<pot.R
        return pot.A*pot.B*exp(-pot.B*r)-6pot.C/r^7
    else
        return 0.0
    end
end

######################################################
# Buckingham Potential                               #
######################################################
