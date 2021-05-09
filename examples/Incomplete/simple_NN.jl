using Pkg; Pkg.activate("."); Pkg.instantiate();

# macro import_these(pkgs...)
#     Pkg.add.(string.(pkgs))
#     text = foldl(*, ["using $p; " for p in pkgs])
#     :($text)
# end
# Pkg.develop(path="../MakeReports")

using Plots, Flux, IterTools, Statistics
using MakeReports
using MDSimulator
report = Report()


text="""
Creating an LJ potential with following paramters:
sigma = 3.4
expilon = 0.99
cutoff = 3*sigma
"""
push!(report, StringPart(text))

# create LJ Potential
σ_ = 3.4
lj_params = LennardJonesParameters(0.99, σ_, 3*σ_)
r = 1.0σ_:0.01σ_:2σ_
𝐄 = [potential_energy(i, lj_params, (1,1)) + 0.05randn() for i in r]

# plot
fig = scatter(r, 𝐄)
savefig(fig, "./output/lj.png")
push!(report, ImagePart("./output/lj.png"; desc="LJ potential."))

text="""
Create data set by normalized energy values.
"""
push!(report, StringPart(text))
# creating normalized dataset
m_ = mean(𝐄)
std_ = std(𝐄)
𝐄 = @. (𝐄 - m_)/std_
new_lj_pot = (x) -> model(x)*std_ +m_
scatter(r, 𝐄)

# creating an MLP model
k = 10
act = σ
model = Chain(x-> [x], Dense(1, k, act), Dense(k, k, act), Dense(k, 1), (x) -> first(x))
ps = Flux.params(model)
ys = model.(r)
text="""
Create an MLP model with following paramters:

activation: sigmoid

hidden layers: [$k, $k, $k]

model: $(model)
"""
push!(report, StringPart(text))

# plot model before training
fig = scatter(r, 𝐄)
plot!(r, ys)
savefig(fig, "./output/before_training.png")
push!(report, ImagePart("./output/before_training.png"; desc="Model before training."))


# set loss function
rx = 2.5σ_:0.1σ_:6σ_
rdata = [(0.0-m_)/std_ for i in rx]

rx2 = 0.75σ_:0.01σ_:σ_
f6 = (x) -> 4*0.99*((σ_/x)^12 - (σ_/x)^6)
rdata2 = [(f6(i)-m_)/std_ for i in rx2]

loss_f(x, y) = sum(@. (model(x)-y)^2) + 0.01*sum(@. ((model(rx2)-rdata2))^2 ) + 1*sum(@. (model(rx)-rdata)^2)

loss_full() = sum([loss_f(x, y) for (x, y) in zip(r, 𝐄)])

# optmizer
opt = ADAM(0.01, (0.9, 0.8))
opt = ADAGrad()
text="""
Using ADAGrad optimizer.
"""
push!(report, StringPart(text))

# training loop
for i in 1:100
    Flux.train!(loss_f, ps, ncycle([(r, 𝐄)], 1000), opt, cb = Flux.throttle(()->(@show loss_full()), 3))
    fig = scatter(r, 𝐄)
    plot!(r, model.(r))
    savefig(fig, "./output/intermediate.png")
end
text="""
Training for 100000 epochs.
"""
push!(report, StringPart(text))


# Final plot
r2 = 0.25σ_:0.01σ_:10σ_
# scatter!(r2, [potential_energy(i, lj_params) for i in r2])
fig = scatter(r2, [potential_energy(i, lj_params, (1,1)) for i in r2], label="LJ")
# scatter!(rx2, rdata2*std_.+m_, color="green")
scatter!(rx, rdata*std_.+m_, color="yellow", label="Zero data")
vline!([0.75σ_], lw=2)
plot!(r2, new_lj_pot.(r2); linewidth=2, label="MLP model")
ylims!(-2, 2)
savefig(fig, "./output/after_training1.png")
ylims!(-2, 100)
savefig(fig, "./output/after_training2.png")

push!(report, ImagePart("./output/after_training1.png"; desc="Model after training (near mimima)."))
push!(report, ImagePart("./output/after_training2.png"; desc="Model after training (repulsive region)."))


# save model
using BSON: @save
@save "./output/trained_NN.bson" model
push!(report, PathPart("./output/trained_NN.bson"; desc="Model after training."))

#print report
print_report("./output/simple_NN", report)
