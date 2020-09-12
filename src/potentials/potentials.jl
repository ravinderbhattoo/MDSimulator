export PotentialParameters, PairPotential, MLPotential, copy_pot

abstract type PotentialParameters end

abstract type PairPotential <: PotentialParameters end
abstract type MLPotential <: PotentialParameters end

include("./lennard_jones.jl")
include("./buckingham.jl")
include("./ML.jl")

abstract type _3BodyPotential <: PotentialParameters end
abstract type _4BodyPotential <: PotentialParameters end
abstract type NBodyPotential <: PotentialParameters end
abstract type BondedPotential <: PotentialParameters end
abstract type FieldPotentials <: PotentialParameters end


function copy_pot(pot)
    args = []
    for f in fieldnames(typeof(pot))
        push!(args, ustrip(getproperty(pot,f)))
    end
    similar(pot, args)
end


copy_pot(pot::MLPotential) = pot
