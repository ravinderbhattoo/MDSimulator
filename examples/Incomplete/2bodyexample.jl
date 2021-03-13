# Imports

using MDSimulator
using Plots
using LinearAlgebra
using ForwardDiff
using Unitful

# Setup
T = 100.0u"K"
m = 39.95u"g/mol"
n = 3
N = n^3
L = uconvert(UNITS.distance, 1.2*n*1u"Å")*3
Δτ = 1.0u"fs"

function rec_grid_(n)
    out = zeros(3, n^3)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                out[:,i + (j-1)*n + (k-1)*n^2] = [i/n, j/n, k/n]
            end
        end
    end
    out
end

u0 = L*rec_grid_(n)
v0 = 0u0/1u"fs"
v0 = maxwell_boltzmann_velocity(N, 100u"K", m)

mass = m*ones(size(u0, 2))
boundary_condition = CubicPBC(ustrip(L))

ϵ = 0.9957u"kJ/mol" # kJ/mol
σpot = 340u"pm" # 3.4Å
R = min(2.5σpot*10, sqrt(3)/2*L)

interatomic_potentials = [LennardJonesParameters([1], [100], [R/1unit(R)])]

u0[1,2] += 0.01u"pm"
# Define simulation object
sim = MDSimulator.SimInfo(u0, v0, mass, interatomic_potentials, boundary_condition, false ; Δτ = Δτ, save_every = 100, thermo_save_every = 100, max_neighs_hard_set=150, reneighboring_every=100)

# Define ensemble
ensemble = [NVE()]

set_thermo(Symbol[:pe, :ke, :temp], ckeys_=String["density"], output="dtab-delim")
set_thermo()

function cal_density(v, x, params)
    return ustrip(uconvert(u"g/cm/cm/cm", (UNITS.mass/CONSTANTS.Avogadro)/UNITS.distance^3 * sum(params.S.sim.mass)/MDBase.volume(params.S.sim.boundary_condition)))
end

# Simulate for n number of steps
sol, parameters = MDSimulator.simulate(10000, sim, ensemble, callbacks=[cb_savethermo_f, cb_reneighboring_f], verbose=true, cthermo=Dict("density"=>cal_density))

# Plots
display(plot_energy(sol, parameters))

display(plot_temperature(sol, parameters))

display(plot_momentum(sol, parameters))






write_trajectory_xyz("./output/ovito/2body", sol, parameters)
write_trajectory_xyz("./output/ovito/2body2", sol2, parameters2)

T_(v) =  0.5*first((sum(v.^2, dims=1)*1ustrip(mass)))

V_(v, u) = MDBase.potential_energy(v, u, [MDSimulator.copy_pot(i) for i in interatomic_potentials], parameters)

∂D_∂v(v, u) = ForwardDiff.gradient((v) -> 𝐷_(v, u), v)

Tᵣ = 1ustrip(3/2*CONSTANTS.kb*1000.101409u"K")
get_temperature(v0, mass)
function 𝐷_(v, u)
    e = 0.0
    m = 1ustrip(mass)
    for i in 1:length(v)
        e += Tᵣ*v[i]^2 /2 - 0.5/4*m[1]*v[i]^4
    end
    return -1.0e-1*e
end

∂D_∂v(1ustrip(v0),1ustrip(u0))


𝐿(v, u) = T_(v) - V_(v, u)

ip = Lagrangian(𝐿)#, 𝐷_)

ip.energy(1ustrip(v0), 1ustrip(u0))
MDBase.get_potential_energy(1ustrip(v0), 1ustrip(u0), parameters) + get_kinetic_energy(1ustrip(v0), parameters.S.sim.mass)

MDBase.get_acceleration(1ustrip(v0), 1ustrip(u0), parameters)

ip.acceleration(1ustrip(v0), 1ustrip(u0))

# Define simulation object
sim2 = MDSimulator.SimInfo(u0, v0, mass, [ip], boundary_condition, false ; Δτ = Δτ, save_every = 10, thermo_save_every = 10, max_neighs_hard_set=150, reneighboring_every=9999999999)

# Define ensemble
ensemble = [NVE()]


set_thermo()
# Simulate for n number of steps
sol2, parameters2 = MDSimulator.simulate(500, sim2, ensemble, callbacks=[cb_savethermo_f, cb_reneighboring_f], verbose=true)

# Plots
display(plot_energy(sol, parameters))
display(plot_energy(sol2, parameters2))

display(plot_temperature(sol, parameters))
display(plot_temperature(sol2, parameters2, 100))

display(plot_momentum(sol, parameters))
display(plot_momentum(sol2, parameters2))


plot([sol(i).x[2][1,:][1] for i in sol.t], sol.t)
plot!([sol(i).x[2][1,:][2] for i in sol.t], sol.t)
plot!([sol(i).x[2][1,:][3] for i in sol.t], sol.t)
plot!([sol(i).x[2][1,:][4] for i in sol.t], sol.t)

plot([sol2(i).x[2][1,:][1] for i in sol2.t], sol2.t)
plot!([sol2(i).x[2][1,:][2] for i in sol2.t], sol2.t)
plot!([sol2(i).x[2][1,:][3] for i in sol2.t], sol2.t)
plot!([sol2(i).x[2][1,:][4] for i in sol2.t], sol2.t)


#

prob, dt, saveat, params = MDSimulator.problem(100, sim, ensemble, callbacks=[cb_savethermo_f, cb_reneighboring_f], verbose=false)

@code_warntype MDBase.solve(prob, MDBase.VelocityVerlet(), dt=sim.Δτ, saveat=saveat)

#
