# export
export potential_energy_plot

# include all potential files here
include("./lennard_jones.jl")
include("./buckingham.jl")
include("./electrostatic.jl")
include("./ML.jl")


"""
    plot_potential_energy_function(pot::T)

Plot potential energy function from PotentialParameter object.

# Arguments:
- `pot`: PotentialParameter object

# Example(s)
```julia
plot_potential_energy_function(LennardJonesParameters())
```
"""
function plot_potential_energy_function(pot::T; init=nothing) where T<:PairPotential
    n = pot.R.s
    fig = Plots.plot()
    for i in 1:n
        for j in i:n
            pair = i,j
            step = pot.R[i,j]/1000
            if init==nothing
                init = step
            end
            r = init:step:pot.R[i,j]
            PE = [potential_energy(i, pot, pair)  for i in r]
            Plots.plot!(1ustrip(r), 1ustrip(PE), label="$i, $j")
        end
    end
    Plots.xlabel!("Distance")
    Plots.ylabel!("Potential energy")
    fig
end


"""
    Create copy of PotentialParameter object after removing attached units.
"""
function copy_pot(pot)
    args = []
    for f in fieldnames(typeof(pot))
        push!(args, 1ustrip(getproperty(pot,f)))
    end
    similar(pot, args)
end

copy_pot(pot::T) where T <: MLPotential = pot
copy_pot(pot::Dummy) = pot

"""
    Unit of energy from potential function.
"""
energy_unit(pot::T) where T = energy_unit(pot)

"""
    Unit of force from potential function.
"""
force_unit(pot::T) where T = force_unit(pot)

"""
    Return potential energy function.
"""
potential_expression(pot::T) where T = potential_expression(pot)
