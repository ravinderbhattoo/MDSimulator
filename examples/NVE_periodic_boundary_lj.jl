using MDSimulator
using Unitful
using ParticlesTools
using ParticlesMesh
using DelimitedFiles
using DiffEqBase
using Statistics

function density2N(ρ, W)
    AvgN = 6.02214e23u"1/mol"
    ρN = uconvert(u"pm^-3", ρ/W*AvgN)
    ρN
end

T = 700.0u"K" # K

ρ = 1449u"g/m^3" # g/m^3
m = 39.95u"g/mol" # g/mol
ρN = density2N(ρ, m) #/A3
L = uconvert(UNITS.distance, 20.0u"Å")
V = L^3
Δτ = 1u"fs"
N = 216
n = Int(round(N^(1/3)))
N = n*n*n
u0 = L*rectangular([1/n, 1/n, 1/n],[n, n, n])
u0 += L/100*randn(size(u0))
v0 = maxwell_boltzmann_velocity(size(u0,2), T, m)
# v0[2:3,:] *= 0

mass = m*ones(size(u0, 2))

uconvert(u"K", get_temperature(v0, mass))

a_ids = ones(Int64,size(u0, 2))

boundary_condition = CubicPBC(ustrip(L))

ϵ = 0.9957u"kJ/mol" # kJ/mol
σ = 340u"pm" # 3.4Å
R = min(2.5σ, 0.5L)

interatomic_potentials = [LennardJonesParameters(ϵ,σ,R)]

Sim = MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids=a_ids, m_ids=a_ids, Δτ = Δτ, save_every = 100, global_save_every = 100, max_neighs_hard_set=150, reneighboring_every=100)

thermostat = AndersenThermostat(300.0)
ensemble = [NVT([thermostat])]
res = simulate(Sim, 5000, ensemble, verbose=true)

Sim2 = sim_init(res, Sim)
thermostat = NoseHooverThermostat(1e6, 216*3, 1, 300.0)
ensemble = [NVT([thermostat])]
res2 = simulate(Sim2, 1000, ensemble, verbose=true)

paraview_save("LJ_system", res)

KE, temp, PE = summary_(res)

TE = KE + PE

using Plots

plot(1:length(KE), [ustrip(i-KE[1]) for i in KE], label="Kinetic Energy")
plot!(1:length(KE), [ustrip(i) for i in PE], label="Potential Energy")
plot!(1:length(KE), [ustrip(i) for i in TE], label="Total Energy")
xlabel!( "Time ($(UNITS.time))")
ylabel!( "Energy ($(UNITS.energy))")


plot(1:length(temp), [ustrip(i) for i in temp], label="Temperature")
xlabel!( "Time ($(UNITS.time))")
ylabel!( "Temperature ($(UNITS.temperature))")

save_trajectory("./output/result1.jld2", res)
res2 = load("./output/result1.jld2")["res"]

res2


#
