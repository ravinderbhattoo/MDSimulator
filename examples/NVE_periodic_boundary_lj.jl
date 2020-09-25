using MDSimulator
using ParticlesMesh

T = 100.0u"K"
m = 39.95u"g/mol"
L = uconvert(UNITS.distance, 15.0u"Å")
Δτ = 1.0u"fs"
N = 125
n = 5

u0 = L*rectangular([1/n, 1/n, 1/n],[n, n, n])
u0 += L/100*randn(size(u0))
v0 = maxwell_boltzmann_velocity(size(u0,2), T, m)

mass = m*ones(size(u0, 2))
boundary_condition = CubicPBC(ustrip(L))

ϵ = 0.9957u"kJ/mol" # kJ/mol
σpot = 340u"pm" # 3.4Å
R = min(2.5σpot, sqrt(3)/2*L)

# interatomic_potentials = [NN(model)]
interatomic_potentials = [LennardJonesParameters(ϵ, σpot, R)]

sim = MDSim(u0, v0, mass, interatomic_potentials, boundary_condition, true; Δτ = Δτ, save_every = 100, thermo_save_every = 100, max_neighs_hard_set=150, reneighboring_every=100)

# thermostat = AndersenThermostat(300.0)
ensemble = [NVE()]


using BenchmarkTools

intg = MDBase.integrator(sim, 400, ensemble, verbose=false)

intg.p.params.M.step = 100

@btime MDSimulator.savethermo_affect_f!(intg)



plot_energy(res)
plot_temperature(res)


a1 = 115528
a2 = 201558
a3 = 287588
a4 = 373618

a4-a3
a3-a2
a2-a1

a1 - (a4-a3)

# Sim2 = sim_init(res, Sim)
# thermostat = NoseHooverThermostat(1e6, 216*3, 1, 300.0)
# ensemble = [NVT([thermostat])]
# res2 = simulate(Sim2, 1000, ensemble, verbose=true)

paraview_save("./output/LJ_system", res)
write_trajectory("./output/ovito/LJ_system", res)

#
