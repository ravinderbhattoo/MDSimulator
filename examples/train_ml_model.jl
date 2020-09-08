using Flux
using JLD2
using DifferentialEquations
using MDSimulator

res = load("./output/result1.jld2")["res"]

function load_data(res; frames=nothing)
    xs = []
    ys = []
    if frames==nothing
        u = res["u"]
        acc = res["acceleration"]
    else
        u = res["u"][frames]
        acc = res["acceleration"][frames]
    end
    for item in u[1:end-1]
        data = Array(reshape(Array(item),(3,:)))'
        N = Int(size(data, 1)/2)
        v = data[1:N,:]
        x = data[N+1:end,:]
        push!(xs, hcat(x, v)')
    end
    for i in 1:size(acc,3)-1
        push!(ys, acc[:,:,i])
    end
    xs, ys
end

xs, ys = load_data(res)

train_loader = Flux.Data.DataLoader((xs, ys), batchsize=2, shuffle=true)


train_loader =  Flux.Data.DataLoader((xs, ys), batchsize=1)

for epoch in 1:100
    for (x, y) in train_loader
        print(x)
        @assert size(x) == (10, 2)
        @assert size(y) == (2,)
    end
end


m = Chain(Dense(6, 15, relu), Dense(15, 14, relu), Dense(14, 3))

ps = Flux.params(m)

loss_f(x, y) = sum((m(x).-y).^2)

loss_full() = sum([loss_f(x, y) for (x,y) in zip(xs,ys)])

opt = ADAM(0.01, (0.9, 0.8))

using IterTools: ncycle
Flux.train!(loss_f, ps, ncycle([(x,y) for (x,y) in zip(xs,ys)], 10), opt, cb = Flux.throttle(()->(@show loss_full()), 3))


#
