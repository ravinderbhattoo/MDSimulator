using MDSimulator
using Unitful
using ParticlesTools
using ParticlesMesh
using DelimitedFiles


function density2N(ρ, W)
    AvgN = 6.02214e23u"1/mol"
    ρN = uconvert(u"Å^-3", ρ/W*AvgN)
    ρN
end

ϵ = 0.9957u"kJ/mol" # kJ/mol
σ = 3.4u"Å" # Å
T = 300.0u"K" # K
kb = 8.3144598e-3u"kJ/K/mol" # kJ/(K*mol)
ϵ = 0.9957u"kJ/mol" # kJ/mol
σ = 3.4u"Å" # Å
ρ = 1449u"g/m^3" # g/m^3
m = 39.95u"g/mol" # g/mol
ρN = density2N(ρ, m) #/A3
L = 100.0u"Å"
V = L^3
N = Int(round(ρN*V))
v_dev = uconvert(u"Å/fs", sqrt(3kb * T / m))
Δτ = 1u"fs"
unit_factor = uconvert(unit(L/Δτ^2), uconvert(u"m/fs^2", 1unit(ϵ/L/m)))


mass = m.val*ones(N)
a_ids = ones(N)

boundary_condition = CubicPBC(L.val)

R = min(6σ, 0.5L)

interatomic_potentials = [LennardJonesParameters(ϵ,σ,R)]
environment = [NVE()]

Sim = MDSim((u0, v0), mass, interatomic_potentials, boundary_condition; a_ids=a_ids, m_ids=a_ids, Δτ = Δτ.val, save_every = 100, global_save_every = 100, max_neighs_hard_set=100, reneighboring_every=100, unit_factor = unit_factor.val)

res = simulate(Sim, 10000, environment)


Mx = WellArray(length(res.t))
Mz = WellArray(length(res.t))
My = WellArray(length(res.t))

for i in 1:length(res.t)
    t = res.t[i]
    data = Array(reshape(Array(res(t)),(3,:)))'
    N = Int(size(data)[1]/2)
    v = data[1:N,:]
    x = data[N+1:end,:]
    M = sum(v, dims=2) / N
    fillit!(Mx, M[1])
    fillit!(My, M[2])
    fillit!(Mz, M[3])
    writedata("./output/sim2_frame_$i.data", (x, v))
end

KE, _, PE = summary_(res)

uconvert(u"kJ/mol", 1unit(m*(L/Δτ)^2))


uf2 = sum([(KE[k+1]-KE[k])/(PE[k+1]-PE[k]) for k in 1:length(KE)-1])/(length(KE)-1)

PE *= -uf2
TE = KE + PE


PE = PE .- (PE[1] - KE[1])
TE = TE .- (TE[1] - KE[1])

using Plots

plot(1:length(KE), [i for i in KE])
plot!(1:length(PE), [i for i in PE])
plot!(1:length(TE), [i for i in TE])


plot(1:length(Mx), [i for i in Mx])
plot!(1:length(My), [i for i in My])
plot!(1:length(Mz), [i for i in Mz])

#
