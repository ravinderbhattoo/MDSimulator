using Flux
using DelimitedFiles
using Statistics

function load_data(folder)
    files = readdir(folder)
    xs = []
    ys = []
    for item in files
        data = readdlm("$folder$(item)", skipstart=2)
        x = data[:,3:5]
        v = data[:,6:8]
        a = data[:,9:11]
        push!(xs, hcat(x, v)')
        push!(ys, a')
    end
    xs, ys
end

folder = "./output/ovito/LJ_system_ovito_files/"
xs, ys = load_data(folder)

train_loader = Flux.Data.DataLoader((xs, ys), batchsize=10, shuffle=true)


train_loader =  Flux.Data.DataLoader((xs, ys), batchsize=1)

function zeroit(m)
    for l in m.layers
        l.W .*= 0.0
        l.b .*= 0.0
    end
end

m1 = Chain(Dense(3, 3, relu), Dense(3, 3, relu))
m2 = Chain(Dense(3, 3, relu), Dense(3, 3, relu))
m3 = Chain(Dense(3, 3, relu), Dense(3, 3))

f(x) = 0.1m1(x[1:3,:]) #+ 0m2(x[4:6,:])

model = Chain(f, m3)

ps = Flux.params(model)

loss_f(x, y) = mean((model(x).-y).^2)

loss_full() = mean([loss_f(x, y) for (x,y) in zip(xs,ys)])

opt = ADAM(0.01, (0.9, 0.8))

# opt = Flux.Optimise.AMSGrad()

using IterTools: ncycle
Flux.train!(loss_f, ps, ncycle([(x,y) for (x,y) in zip(xs,ys)], 100), opt, cb = Flux.throttle(()->(@show loss_full()), 3))




#
