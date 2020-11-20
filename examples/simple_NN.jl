using Plots
using Flux
using IterTools: ncycle
using Statistics
using MDSimulator

σ_ = 3.4
lj_params = LennardJonesParameters(0.99, σ_, 3*σ_)
r = 1.0σ_:0.01σ_:2σ_
𝐄 = [potential_energy(i, lj_params) + 0.05randn() for i in r]
scatter(r, 𝐄)
# 𝐄 = [f(i) for i in r]

m_ = mean(𝐄)
std_ = std(𝐄)
𝐄 = @. (𝐄 - m_)/std_

new_lj_pot = (x) -> model(x)*std_ +m_

scatter(r, 𝐄)

k = 10
act = σ
model = Chain(x-> [x], Dense(1, k, act), Dense(k, k, act), Dense(k, 1), (x) -> first(x))
ps = Flux.params(model)
ys = model.(r)

scatter(r, 𝐄)
plot!(r, ys)

rx = 2.5σ_:0.1σ_:6σ_
rdata = [(0.0-m_)/std_ for i in rx]

rx2 = 0.75σ_:0.01σ_:σ_
f6 = (x) -> 4*0.99*((σ_/x)^12 - (σ_/x)^6)
rdata2 = [(f6(i)-m_)/std_ for i in rx2]

loss_f(x, y) = sum(@. (model(x)-y)^2) + 0.01*sum(@. ((model(rx2)-rdata2))^2 ) + 1*sum(@. (model(rx)-rdata)^2)

loss_full() = sum([loss_f(x, y) for (x, y) in zip(r, 𝐄)])

opt = ADAM(0.01, (0.9, 0.8))
opt = ADAGrad()

Flux.train!(loss_f, ps, ncycle([(r, 𝐄)], 100000), opt, cb = Flux.throttle(()->(@show loss_full()), 3))

scatter(r, 𝐄)
plot!(r, model.(r))


r2 = 0.75σ_:0.01σ_:10σ_
plot(r2, new_lj_pot.(r2))
# scatter!(r2, [potential_energy(i, lj_params) for i in r2])
scatter!(r2, [potential_energy(i, lj_params) for i in r2], color="red")
# scatter!(rx2, rdata2*std_.+m_, color="green")
scatter!(rx, rdata*std_.+m_, color="yellow")
vline!([0.75σ_], lw=2)
#

using BSON: @save
@save "./output/trained_NN.bson" model
