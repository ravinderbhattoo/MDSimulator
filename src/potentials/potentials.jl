# include all potential files here

include("./lennard_jones.jl")
include("./buckingham.jl")
include("./ML.jl")

function copy_pot(pot)
    args = []
    for f in fieldnames(typeof(pot))
        push!(args, 1ustrip(getproperty(pot,f)))
    end
    similar(pot, args)
end

copy_pot(pot::T) where T <: Union{MLPotential,Zygote.Params} = pot
copy_pot(pot::Dummy) = pot
