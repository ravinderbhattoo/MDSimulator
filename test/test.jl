using MDSimulator
using Unitful
using ParticlesTools
using ParticlesMesh
using DelimitedFiles
using DiffEqBase

function density2N(ρ, W)
    AvgN = 6.02214e23u"1/mol"
    ρN = uconvert(u"pm^-3", ρ/W*AvgN)
    ρN
end

T = 300.0u"K" # K

ρ = 1449u"g/m^3" # g/m^3
m = 39.95u"g/mol" # g/mol
ρN = density2N(ρ, m) #/A3
L = uconvert(UNITS.distance, 500.0u"Å")
V = L^3
N = Int(round(ρN*V))

Δτ = 1u"fs"
N1, N2, N3 = 13, 14, 15
u0 = L*rectangular([1/N3,1/N2,1/N3],[N1, N2, N3])
u0 = L/100*randn(size(u0)) .+ u0

v0 = maxwell_boltzmann_velocity(N, T, m)

mass = m*ones(N)
a_ids = ones(Int64,N)

boundary_condition = CubicPBC(ustrip(L))

ϵ = 0.9957u"kJ/mol" # kJ/mol
σ = 340u"pm" # 3.4Å
R = min(6σ, 0.5L)

interatomic_potentials = [LennardJonesParameters(ϵ,σ,R)]
environment = [NVE()]

Sim = MDSim(u0, v0, mass, interatomic_potentials, boundary_condition; a_ids=a_ids, m_ids=a_ids, Δτ = Δτ, save_every = 1, global_save_every = 100, max_neighs_hard_set=100, reneighboring_every=100)

res = simulate(Sim, 10000, environment, verbose=true)

KE, temp, PE = summary_(res)

TE = KE + PE

using Plots
plot(res.t[2:end], [ustrip(i-KE[1]) for i in KE], label="Kinetic Energy")
plot!(res.t[2:end], [ustrip(i-PE[1]) for i in PE], label="Potential Energy")
plot!(res.t[2:end], [ustrip(i-TE[1]) for i in TE], label="Total Energy")
xlabel!( "Time ($(UNITS.time))")
ylabel!( "Energy ($(UNITS.energy))")


plot(res.t[2:end], [ustrip(i) for i in temp], label="Temperature")
xlabel!( "Time ($(UNITS.time))")
ylabel!( "Temperature ($(UNITS.temperature))")

save_trajectory("/tmp/result.jld", res)
res2 = load("/tmp/result.jld")["res"]

#
