using ForwardDiff
using LinearAlgebra

# exports
export potential_force, potential_energy

######################################################
# NN Potential                                       #
######################################################

abstract type MLPairPotential <: MLPotential end
abstract type MLModelPotential <: MLPotential end
export MLPairPotential

export NNPair
struct NNPair <: MLPairPotential
    model
    force
    energy
    function NNPair(model)
        function energy(x)
            return model(x)
        end
        neg_∂U_∂x = (x) -> -ForwardDiff.derivative(energy, x)
        new(model, neg_∂U_∂x, energy)
    end
end

function Base.show(stream::IO, pp::NNPair)
    println(stream, "NN Pair Potential:")
    println(stream, "\tPotential Energy (𝑈) = model(x) " )
    println(stream, "\tPotential force = -∂𝑈/∂x where" )
    print(stream, "\tmodel:\t"); show(stream, pp.model); println(stream)
end


function potential_energy(r, pot::NNPair) 
    return pot.energy(r)
end

function potential_force(r, pot::NNPair)
    return pot.force(r)
end

######################################################
# NN Potential                                       #
######################################################
######################################################
# LNN Potential                                       #
######################################################

export LNN 

struct LNN <: MLModelPotential
    model
    acceleration
    energy
    function LNN(model)
        𝐿(ẋ, x) = reshape([sum(model(vcat(x,ẋ)))], (1,1))
        ∇ₓ(ẋ, x, 𝐿) = ForwardDiff.jacobian((x) -> 𝐿(ẋ, x), x)
        ∇ₓ̇(ẋ, x, 𝐿) = ForwardDiff.jacobian((ẋ) -> 𝐿(ẋ, x), ẋ)
        ∇ₓ̇ᵀ(ẋ, x, 𝐿) = transpose(∇ₓ̇(ẋ, x, 𝐿))
        ∇ₓ̇∇ₓ̇ᵀ(ẋ, x, 𝐿) = ∇ₓ̇(ẋ, x, (ẋ, x) -> ∇ₓ̇ᵀ(ẋ, x, 𝐿))
        ∇ₓ∇ₓ̇ᵀ(ẋ, x, 𝐿) = ∇ₓ(ẋ, x, (ẋ, x) -> ∇ₓ̇ᵀ(ẋ, x, 𝐿))
        function acc(ẋ, x)
            ∇ₓ𝐿 = ∇ₓ(ẋ, x, 𝐿)'
            ∇ₓ̇∇ₓ̇ᵀ𝐿 = ∇ₓ̇∇ₓ̇ᵀ(ẋ, x, 𝐿)
            ∇ₓ∇ₓ̇ᵀ𝐿 = ∇ₓ∇ₓ̇ᵀ(ẋ, x, 𝐿)
            N = length(x)
            q̈ = reshape(inv(∇ₓ̇∇ₓ̇ᵀ𝐿 + 1e-9I(N))*(∇ₓ𝐿 - ∇ₓ∇ₓ̇ᵀ𝐿 * reshape(ẋ, (:,1))), size(x)...)
            return q̈
        end
        new(model, acc, (ẋ, x)->first(𝐿(ẋ, x)))
    end
end


function Base.show(stream::IO, pp::LNN)
    println(stream, "LNN Potential:")
    println(stream, "\tPotential Energy = model(x) where" )
    print(stream, "\tmodel:\t"); show(stream, pp.model); println(stream)
end


function potential_energy(v::Array{F1,2}, u::Array{F2,2}, pot::LNN) where {F1, F2} 
    return pot.energy(v, u)
end

function acceleration(v::Array{F1,2}, u::Array{F2,2}, pot::LNN) where {F1, F2}
    return pot.acceleration(v, u)
end

######################################################
# LNN Potential                                       #
######################################################
