# exports
export potential_force, potential_energy, PotentialParameters

######################################################
# Buckingham Potential                               #
######################################################

export BuckinghamParameters

struct BuckinghamParameters{AType <: Number, BType <: Number, CType <: Number, lType <: Number} <: PairPotential
    A::AType
    B::BType
    C::CType
    R::lType
    equillibrium::lType
    well_energy::AType
    energy_shift::AType
end

function BuckinghamParameters(A::Number, B::Number, C::Number, R::Number)
    energy_shift = A*exp(-B*R) - C/R^6
    r = R
    for i in 1:100
        r = (6C/(A*B)*exp(-B*r))^(1/7)
    end
    well_energy = A*exp(-B*r) - C/r^6
    BuckinghamParameters(A, B, C, R, r, well_energy, -energy_shift)
end

function BuckinghamParameters()
    BuckinghamParameters(1.0, 1.0, 0.001, 2.5)
end

function Base.show(stream::IO, pp::BuckinghamParameters)
    println(stream, "Buckingham:")
    println(stream, "\tV(r) = A e⁻ᴮʳ - C r⁻⁶ where")
    print(stream, "\tA:\t"); show(stream, pp.A); println(stream)
    print(stream, "\tB:\t"); show(stream, pp.B); println(stream)
    print(stream, "\tC:\t"); show(stream, pp.C); println(stream)
    print(stream, "\tR:\t"); show(stream, pp.R); println(stream)
    print(stream, "\tEquillibrium position:\t"); show(stream, pp.equillibrium); println(stream)
    print(stream, "\tWell energy:\t"); show(stream, pp.well_energy); println(stream)
    print(stream, "\tEnergy shift:\t"); show(stream, pp.energy_shift); println(stream)
end



function potential_energy(r::Number, pot::BuckinghamParameters)
    if r<pot.R
        return pot.A*exp(-pot.B*r) - pot.C/r^6 + pot.energy_shift
    else
        return 0.0
    end
end

function potential_force(r::Number, pot::BuckinghamParameters)
    if r<pot.R
        return pot.A*pot.B*exp(-pot.B*r)-6pot.C/r^7
    else
        return 0.0
    end
end

######################################################
# Buckingham Potential                               #
######################################################


#
