using MDSimulator
using ParticlesTools
using ParticlesMesh
using DelimitedFiles


u0 = 2.0*rectangular([1.0,1,1],[10,10,10])
v0 = 0.1*randn(size(u0))
v0 .-= sum(v0, dims=2)/size(v0)[2]

mass = ones(size(u0)[2])
a_ids = ones(size(u0)[2])
boundary_condition = CubicPBC(20.0)

interatomic_potentials = [LennardJonesParameters(1.0e-3,1.0,5.0)]
environment = [NVE()]

Sim = MDSim((u0, v0), mass, interatomic_potentials, boundary_condition; a_ids=a_ids, m_ids=a_ids, Δτ = 1.0e-3, save_every = 100, global_save_every = 10, max_neighs_hard_set=200, reneighboring_every=100)

res = simulate(Sim, 1000, environment)


for i in 1:length(res.t)
    t = res.t[i]
    data = Array(reshape(Array(res(t)),(3,:)))'
    N = Int(size(data)[1]/2)
    v = data[1:N,:]
    x = data[N+1:end,:]
    writedata("./output/frame_$i.data", (x, v))
end

KE, _ = summary(res)

using Plots

plot(1:length(KE), [i for i in KE])


#
