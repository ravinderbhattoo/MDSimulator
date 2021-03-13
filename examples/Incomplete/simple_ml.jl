using Flux
using Flux.Optimise: update!
using Plots

xs = 400:10:500

struct lj
    σ::Array{Float64, 1}
    ϵ::Array{Float64, 1}
end

energy(r, pot) = first(@. 4pot.ϵ*((pot.σ/r)^12 - (pot.σ/r)^6))
energy(r, pot::Chain) = first(pot([r]))

LJ = lj([340], [1])
data = [energy(x, LJ) for x ∈ xs]

LJ2 = lj([400], [2])
Flux.trainable(a::lj) = (a.σ, a.ϵ)

ps = params(LJ2)

loss(xs, data) = sum(([energy(x, LJ2) for x ∈ xs] - data).^2)

opt = ADAGrad(10)
using IterTools: ncycle
Flux.train!(loss, ps, ncycle([(xs, data)], 100), opt, cb = make_plot)


function make_plot()
  println("\tLoss: ", loss())
  ys = [energy(x, LJ2) for x ∈ xs]
  data1 = [energy(x, LJ) for x ∈ xs]
  fig = scatter(xs, data1, label="Data")
  plot!(xs, ys, label="Model")
  display(fig)
end


make_plot()

# #
