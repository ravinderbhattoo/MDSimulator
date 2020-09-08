abstract type PotentialParameters end

abstract type PairPotential <: PotentialParameters end

include("./lennard_jones.jl")
include("./buckingham.jl")

abstract type _3BodyPotential <: PotentialParameters end
abstract type _4BodyPotential <: PotentialParameters end
abstract type NBodyPotential <: PotentialParameters end
abstract type BondedPotential <: PotentialParameters end
abstract type FieldPotentials <: PotentialParameters end
