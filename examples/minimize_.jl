# imports
using MDSimulator
using ParticlesMesh

# Setup
T = 100.0u"K"
m = 39.95u"g/mol"
L = uconvert(UNITS.distance, 10u"Å")
Δτ = 1.0u"fs"
N = 1000
n = 3

u0 = L*rectangular([1/n, 1/n, 1/n],[n, n, n])
u0 += L/100*randn(size(u0))
v0 = maxwell_boltzmann_velocity(size(u0,2), T, m)

mass = m*ones(size(u0, 2))
boundary_condition = CubicPBC(ustrip(L))

ϵ = 0.9957u"kJ/mol" # kJ/mol
σpot = 340u"pm" # 3.4Å
R = min(2.5σpot, sqrt(3)/2*L)

interatomic_potentials = [LennardJonesParameters(ϵ, σpot, R)]

# Define simulation object
sim = MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, true; Δτ = Δτ, save_every = 1000, thermo_save_every = 1000, max_neighs_hard_set=150, reneighboring_every=100,)

# Define ensemble for simulation
ensemble = [NVE()]

parameters = MDSimulator.exe_at_start(sim, 1, true)
println("Minimum PE: ", get_potential_energy(sim.v0, sim.u0, parameters, sim))

# Minimize energy
sim = minimize(sim)

parameters = MDSimulator.exe_at_start(sim, 1, true)
print("Minimum PE: ", get_potential_energy(sim.v0,sim.u0, parameters, sim))
