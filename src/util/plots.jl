using Plots
using Statistics

#exports
export plot_temperature, plot_energy, plot_momentum

function plot_energy_(res::MDBase.ODESolution; E="PE")
    KE, temp, PE = summary_(res)
    TE = KE + PE
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]

    if "KE"==E
        plot(t, [ustrip(i) for i in KE], label="Kinetic Energy")
    elseif "PE"==E
        plot(t, [ustrip(i) for i in PE], label="Potential Energy")
    elseif "TE"==E
        plot(t, [ustrip(i) for i in TE], label="Total Energy")
    else
        throw("Only PE, KE and TE are allowed.")
    end

    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Energy ($(UNITS.energy))")
end

function plot_energy(res::MDBase.ODESolution, params; E=["PE", "KE", "TE"], steps=0)
    KE, temp, PE = summary_(params)
    TE = KE + PE
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]
    if steps!=0 && steps<N
        KE = KE[1:steps]
        TE = TE[1:steps]
        PE = PE[1:steps]
        t = t[1:steps]
    end
    fig = plot()
    if "KE" in E
        plot!(t, [ustrip(i) for i in KE], label="Kinetic Energy")
    end
    if "PE" in E
        plot!(t, [ustrip(i) for i in PE], label="Potential Energy")
    end
    if "TE" in E
        plot!(t, [ustrip(i) for i in TE], label="Total Energy")
    end
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Energy ($(UNITS.energy))")
    fig
end

function plot_temperature(res::MDBase.ODESolution, params)
    KE, temp, PE = summary_(params)
    N = length(KE)
    fig = plot(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(i) for i in temp], label="Temperature")
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Temperature ($(UNITS.temperature))")
    fig
end

function plot_temperature(res::MDBase.ODESolution, params, window)
    rw = convert(Int64, floor(window/2))
    KE, temp, PE = summary_(params)
    N = length(KE)
    fig = plot(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(mean(temp[max(i-rw,1):min(i+rw,N)])) for i in 1:N], label="Temperature")
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Temperature ($(UNITS.temperature))")
    fig
end


function plot_momentum(res::MDBase.ODESolution, params)
    p = get_momentum_all(res, params)
    N = length(res.t)
    xs = res.t[1]:res.t[end]/(N-1):res.t[end]
    fig = plot(xs, [ustrip(i) for i in p[1,:]], label="P₁")
    plot!(xs, [ustrip(i) for i in p[2,:]], label="P₂")
    plot!(xs, [ustrip(i) for i in p[3,:]], label="P₃")
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Momentum ($(UNITS.velocity*UNITS.mass))")
    fig
end
