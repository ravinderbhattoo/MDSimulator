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
        ∇ₓ!(out, ẋ, x) = ForwardDiff.gradient!(out, (x) -> 𝐿(ẋ, x), x)
        ∇ₓ̇!(out, ẋ, x) = ForwardDiff.gradient!(out, (ẋ) -> 𝐿(ẋ, x), ẋ)
        ∇ₓ̇(ẋ, x) = ForwardDiff.gradient((ẋ) -> 𝐿(ẋ, x), ẋ)
        ∇ₓ̇𝐷!(out, ẋ, x) = ForwardDiff.gradient!(out, (ẋ) -> 𝐷(ẋ, x), ẋ)
        ∇ₓ̇ᵀ!(out, ẋ, x) = transpose(∇ₓ̇!(out, ẋ, x))
        ∇ₓ̇ᵀ(ẋ, x) = ∇ₓ̇(ẋ, x)
        ∇ₓ̇∇ₓ̇ᵀ!(out, ẋ, x) = ForwardDiff.jacobian!(out, (ẋ) -> ∇ₓ̇ᵀ(ẋ, x), ẋ)
        ∇ₓ∇ₓ̇ᵀ!(out, ẋ, x) = ForwardDiff.jacobian!(out, (x) -> ∇ₓ̇ᵀ(ẋ, x), x)
        function acc(ẋ, x)
            N = size(x, 2)
            q̈ = similar(x)
            A_ = rand(3N, 3N)
            B_ = similar(x)
            C_ = rand(3N, 3N)
            D_ = similar(x)
            ∇ₓ̇∇ₓ̇ᵀ!(A_, ẋ, x)
            ∇ₓ!(B_, ẋ, x)
            ∇ₓ∇ₓ̇ᵀ!(C_, ẋ, x)
            ∇ₓ̇𝐷!(D_, ẋ, x)
            function f(A_, B_, C_, D_, N, ẋ, q̈)
                A(i) = @view A_[3i-2:3i,3i-2:3i]
                C(i) = @view C_[3i-2:3i,3i-2:3i]
                B(i) = @view B_[:,i]
                D(i) = @view D_[:,i]
                E(i) = @view ẋ[:,i]
                for i in 1:N
                    q̈[:,i] .= inv(A(i))*(B(i) - D(i) - C(i)*E(i))
                end
                return q̈
            end
            f(A_, B_, C_, D_, N, ẋ, q̈)
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
