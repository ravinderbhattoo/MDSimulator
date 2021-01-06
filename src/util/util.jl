export cutoff, reneighboring!
export cal_pdf, cal_pdfs, plot_pdf, plot_pdf!


cutoff(pot::Any) = maximum(pot.R)
cutoff(pot::MLPotential) = 0.0UNITS.distance
cutoff(pot::Dummy) = 0.0UNITS.distance

function reneighboring!(x, params)
    neigh_mat = find_neighbors(x, params.S.others.cut_dist, params.S.sim.boundary_condition, params.S.sim.others.max_neighs, hard_max = params.S.sim.others.max_neighs_hard_set)
    rows = size(neigh_mat, 1)
    fill!(params.S.others.neighs, 0)
    nrows = min(rows, size(params.S.others.neighs, 1))
    if nrows < rows
        print("Neighs has exceeded($rows) capacity($nrows)\n")
    end
    params.S.others.neighs[1:nrows,:] = neigh_mat[1:nrows,:]
end

function cal_pdf(x, bc, neigh, r2; cutoff=100.0, pair=nothing, a_ids=Int64[], bins=1000)
    r = []
    Nα, Nβ = 0, 0
    if pair != nothing
        if length(a_ids) != 0
            α,β = pair
            Nα = sum(a_ids .== α)
            Nβ = sum(a_ids .== β)
            for βj in 1:size(neigh, 2)
                for αi in 2:size(neigh, 1)
                    if a_ids[βj]==β
                        kk = neigh[αi, βj]
                        if kk!=0
                            if a_ids[kk]==α
                                push!(r, sqrt(r2[αi, βj]))
                            end
                        end
                    end
                end
            end
        else
            @error "Please provide atom ids."
        end
    else
        r = sqrt.(r2)
        Nα = size(x,2)
        Nβ = size(x,2)
    end
    r = [i for i in r if i!=0.0]
    V = MDBase.volume(bc)
    ρα = Nα/V
    dr = cutoff/bins
    h = zeros(Float64, bins)
    for i in r
        ind = convert(Int64, ceil(i/dr))
        h[ind] += (1/ρα)*(1/(dr*4pi*(ind*dr-dr/2)^2))
    end
    x = dr/2:dr:cutoff-dr/2
    h/Nβ, x
end

function cal_pdfs(X, bc, args...; pairs=nothing, cutoff=100.0, kwargs...)
    H = [[0.0] for _ in pairs]
    x_ = [[0.0] for _ in pairs]
    for x in X
        neigh, r2 = MDBase.find_neighbors(x, convert(Float64, cutoff), bc, max(size(x)...), with_distance2=true)
        for j in 1:length(pairs)
            h, x__ = cal_pdf(x, bc, neigh, r2, args...; pair=pairs[j], cutoff=cutoff, kwargs...)
            h .+= H[j]
            H[j] = h
            x_[j] = x__
        end
    end
    H/length(X), x_
end


function plot_pdf!(X, bc, args...; cutoff::Float64=100.0, pcutoff=0.0, plot_prop=Dict(), kwargs...)
    H, x = cal_pdfs(X, bc, args...; cutoff=cutoff, kwargs...)
    if pcutoff==0.0
        pcutoff=cutoff
    end

    expand_(a, i) = if length(a)==i a else collect((a[1] for j in 1:i)) end

    plot_prop = Dict([k=>expand_(plot_prop[k], length(H)) for k in keys(plot_prop)]...)
    for i in 1:length(H)
        Plots.plot!(x[i], H[i]; [k=>plot_prop[k][i] for k in keys(plot_prop)]...)
    end
    Plots.xlims!(0.0, pcutoff)
end


function plot_pdf(X, bc, args...; cutoff::Float64=100.0, pcutoff=0.0, plot_prop=Dict(), kwargs...)
    fig = Plots.plot()
    plot_pdf!(X, bc, args...; cutoff=cutoff, pcutoff=pcutoff, plot_prop=plot_prop, kwargs...)
end
