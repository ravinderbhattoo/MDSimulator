# Imports

using MDSimulator
using ParticlesMesh
using Plots

# Setup
T = 100.0u"K"
m = 39.95u"g/mol"
L = uconvert(UNITS.distance, 15u"Å")
Δτ = 1.0u"fs"

u0 = L*rectangular([1/3, 1, 1],[3, 1, 1])
v0 = 0u0/1u"fs"

mass = m*ones(size(u0, 2))
boundary_condition = CubicPBC(1000000.0)

ϵ = 0.9957u"kJ/mol" # kJ/mol
σpot = 340u"pm" # 3.4Å
R = min(2.5σpot, sqrt(3)/2*L)

interatomic_potentials = [LennardJonesParameters([ϵ, ϵ], [σpot, 0.8σpot], [R, R])]


# Define simulation object
sim = MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, true ; a_ids=[1, 2, 1], Δτ = Δτ, save_every = 10, thermo_save_every = 10, max_neighs_hard_set=150, reneighboring_every=100)

# Define ensemble
ensemble = [NVE()]


# Simulate for n number of steps
sol, parameters = MDSimulator.simulate(50000, sim, ensemble, callbacks=[cb_savethermo_f, cb_reneighboring_f], verbose=true)

# Plots
display(plot_energy(sol, parameters))
display(plot_temperature(sol, parameters))
display(plot_momentum(sol, parameters))


plot([sol(i).x[2][1,:][1] for i in sol.t], sol.t)
plot!([sol(i).x[2][1,:][2] for i in sol.t], sol.t)
plot!([sol(i).x[2][1,:][3] for i in sol.t], sol.t)


#write_trajectory_xyz("./output/ovito/2body", sol, parameters)