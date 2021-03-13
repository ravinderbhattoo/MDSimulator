macro LJ(σ, ϵ)
    return :( (r) -> begin
        temp = ($(esc(σ))./r)^6
        return 4.0 .*$(esc(ϵ)) .* ( temp^2 .- temp )
        end)
end


function lj_s(r, σ, ϵ, ind)
    temp = (σ[ind]./r)^6
    return 4*ϵ[ind] .* ( temp^2 .- temp )
end

lj = @LJ 1 1

r = rand(3, 100)
t = Int64.(rand(1, 100) .< 0.5)

σ = [1]
ϵ = [1]
function PE(r, t, func, σ, ϵ)
    pe = 0.0
    for i in 1:size(r, 2)
        for j in i+1:size(r, 2)
            d = sum((r[:,i] .- r[:,j]).^2)
            if t[i]==1
                pe += func(d, σ, ϵ, t[i])
            else
                pe += 0.0
            end
        end
    end
    return pe
end

using BenchmarkTools

@btime PE(r, t, lj_s, σ, ϵ)

function PE2(r, t, funcs)
    pe = 0.0
    for i in 1:size(r, 2)
        for j in i+1:size(r, 2)
            d = sum((r[:,i] .- r[:,j]).^2)
            if t[i]==1
                pe += funcs(d)
            else
                pe += 0.0
            end
        end
    end
    return pe
end

lj1 = @LJ 1 1
lj2 = @LJ 1 1

@btime PE2(r, t, lj1)

dump(:(1+1))

#
