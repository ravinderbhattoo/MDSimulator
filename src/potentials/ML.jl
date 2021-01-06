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
    function NNPair(model)
        neg_∂U_∂x(m) = (x) -> -ForwardDiff.derivative(m(x), x)
        new(model, SymMat(neg_∂U_∂x.(model)))
    end
end

function Base.show(stream::IO, pp::NNPair)
    println(stream, "NN Pair Potential:")
    println(stream, "\tPotential Energy (𝑈) = model(x) " )
    println(stream, "\tPotential force = -∂𝑈/∂x where" )
    print(stream, "\tmodel:\t"); show(stream, pp.model); println(stream)
end


function potential_energy(r, pot::NNPair, pair)
    return pot.model[pair...](r)
end

function potential_force(r, pot::NNPair, pair)
    return pot.force[pair...](r)
end

######################################################
# NN Potential                                       #
######################################################
######################################################
# LNN Potential                                       #
######################################################

export Lagrangian

struct Lagrangian <: MLModelPotential
    𝐿
    acceleration
    energy
    function Lagrangian(𝐿, 𝐷)
        ∇ₓ(ẋ, x) = ForwardDiff.gradient((x) -> 𝐿(ẋ, x), x)
        ∇ₓ̇(ẋ, x) = ForwardDiff.gradient((ẋ) -> 𝐿(ẋ, x), ẋ)
        ∇ₓ̇𝐷(ẋ, x) = ForwardDiff.gradient((ẋ) -> 𝐷(ẋ, x), ẋ)
        ∇ₓ̇ᵀ(ẋ, x) = transpose(∇ₓ̇(ẋ, x))
        ∇ₓ̇∇ₓ̇ᵀ(ẋ, x) = ForwardDiff.jacobian((ẋ) -> ∇ₓ̇ᵀ(ẋ, x), ẋ)
        ∇ₓ∇ₓ̇ᵀ(ẋ, x) = ForwardDiff.jacobian((x) -> ∇ₓ̇ᵀ(ẋ, x), x)
        function acc(ẋ, x)
            q̈ = 0*x
            N = size(x, 2)
            A_ = ForwardDiff.jacobian((ẋ) -> ∇ₓ̇(ẋ, x), ẋ)
            C_ = ForwardDiff.jacobian((x) -> ∇ₓ̇(ẋ, x), x)
            B_ = ∇ₓ(ẋ, x)
            D_ = ∇ₓ̇𝐷(ẋ, x)
            A(i) = A_[3i-2:3i,3i-2:3i]
            C(i) = C_[3i-2:3i,3i-2:3i]
            B(i) = B_[:,i]
            D(i) = D_[:,i]
            for i in 1:N
                q̈[:,i] .= inv(A(i))*(B(i) - D(i) - C(i)*ẋ[:,i])
            end
            return q̈
        end
        new(𝐿, acc, 𝐿)
    end
end

function Lagrangian(𝐿)
    𝐷(ẋ, x) = 0.0
    Lagrangian(𝐿, 𝐷)
end


function Base.show(stream::IO, pp::Lagrangian)
    println(stream, "LNN Potential:")
    println(stream, "\tPotential Energy = model(x) where" )
    print(stream, "\tmodel:\t"); show(stream, pp.𝐿); println(stream)
end


function potential_energy(v::Array{F1,2}, u::Array{F2,2}, pot::Lagrangian) where {F1, F2}
    return -pot.energy(0v, u)
end

function acceleration(v::Array{F1,2}, u::Array{F2,2}, pot::Lagrangian) where {F1, F2}
    return pot.acceleration(v, u)
end

######################################################
# LNN Potential                                       #
######################################################
