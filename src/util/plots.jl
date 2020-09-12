using Plots

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

function plot_energy(res::MDBase.ODESolution; E=["PE", "KE", "TE"])
    KE, temp, PE = summary_(res)
    TE = KE + PE
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]
    plot()
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
end

function plot_temperature(res::MDBase.ODESolution)
    KE, temp, PE = summary_(res)
    N = length(KE)
    plot(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(i) for i in temp], label="Temperature")
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Temperature ($(UNITS.temperature))")
end


function plot_momentum(res::MDBase.ODESolution)
    p = get_momentum_all(res)
    N = length(res.t)
    plot(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(i) for i in p[1,:]], label="P₁")
    plot!(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(i) for i in p[2,:]], label="P₂")
    plot!(res.t[1]:res.t[end]/(N-1):res.t[end], [ustrip(i) for i in p[3,:]], label="P₃")
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Momentum ($(UNITS.velocity*UNITS.mass))")
end
