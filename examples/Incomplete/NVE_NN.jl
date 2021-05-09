# Imports
using DiffEqFlux
using IterTools: ncycle
using MDSimulator
using Plots
using Flux
using LinearAlgebra
using ForwardDiff
using Unitful
using DiffEqSensitivity
using PDMesh


# Setup
T = 100.0u"K"
m = 39.95u"g/mol"
L = uconvert(UNITS.distance, 15u"Å")
Δτ = 1.0u"fs"

u0 = L/2 *create(Cuboid([0 2;0 2;0 2]), resolution=1.0, rand_=0.05)[1]
v0 = 0u0/1u"fs"

N = size(u0, 2)

mass = m*ones(size(u0, 2))
open_bounadry = OpenBoundary()
boundary_condition = GeneralBC(open_bounadry, open_bounadry, open_bounadry)

ϵ = 0.9957u"kJ/mol" # kJ/mol
σpot = 340u"pm" # 3.4Å
R = min(2.5σpot, sqrt(3)/2*L)

interatomic_potentials = [LennardJonesParameters(ϵ, σpot, R)]

# Define simulation object
sim = MDSimulator.SimInfo(u0, v0, mass, interatomic_potentials, boundary_condition, true; Δτ = Δτ, save_every = 100, thermo_save_every = 100, max_neighs_hard_set=150, reneighboring_every=1_000_000,)

# Define ensemble
ensemble = [NVE()]

# Callbacks
callbacks = [cb_reneighboring_f, cb_savethermo_f]

# Simulate for n number of steps
sol, parameters = MDSimulator.simulate(10000, sim, ensemble, callbacks=callbacks, verbose=true)
   
ẍdata = 1parameters.S.acc
xdata = [1i.x[2] for i in sol.u]



model = Chain((x)->[x], Dense(1, 3, σ), Dense(3, 1, σ), (x)->first(x))
pot1 =  NNPair([model])

# define function to get cutoff
MDSimulator.cutoff(pot::NNPair) = R

θ, re = Flux.destructure(pot1.energy[1,1])

sim2 = MDSimulator.SimInfo(u0, v0, mass, [pot1], boundary_condition, false; Δτ = Δτ, save_every = 100, thermo_save_every = 100, max_neighs_hard_set=150, reneighboring_every=1_000_000,)
# parameters2 = MDSimulator.exe_at_start(sim2, 10000, false)

sim2.interatomic_potentials[1][1] = [pot1]
sol2, parameters2 = MDSimulator.simulate(10000, sim2, ensemble, callbacks=callbacks, verbose=false)

function loss(θ)
    model = re(θ)
    sim2.interatomic_potentials[1][1].energy[1,1].layers[2].W .= model.layers[2].W
    sim2.interatomic_potentials[1][1].energy[1,1].layers[2].b .= model.layers[2].b
    sim2.interatomic_potentials[1][1].energy[1,1].layers[3].W .= model.layers[3].W
    sim2.interatomic_potentials[1][1].energy[1,1].layers[3].b .= model.layers[3].b
    sol2, parameters2 = MDSimulator.simulate(10000, sim2, ensemble, callbacks=callbacks, verbose=false)
    l = 0.0
    for i in length(sol2.u)
        l += sum((sol.u[i].x[2] .- sol2.u[i].x[2]).^2)
    end
    l
end

ϵ₁ = 0.1
η = 0.0001
for i in 1:100
    new_θ = ones(length(θ))*θ' + I(length(θ))*ϵ₁
    l = loss(θ)
    gs = [min(max((loss(new_θ[i,:]) - l)/ϵ₁*η, -0.01*θ[i]), 0.01*θ[i]) for i in 1:length(θ)]
    θ .-= gs
    @show i l
end


r = 0.5:100:3000

[MDSimulator.potential_force(i, sim.interatomic_potentials[1][1]) for i in r]
[MDSimulator.potential_force(i, sim2.interatomic_potentials[1][1]) for i in r]


#

m = ones(100)
v = 2ones(3,100)

for i in v
    println(i)
end

T = 0.5m*v2






#