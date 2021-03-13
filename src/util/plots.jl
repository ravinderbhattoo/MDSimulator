using Plots
using Statistics

#exports
export plot_temperature, plot_energy, plot_momentum

function plot_energy(res::MDBase.ODESolution, params; E=["PE", "KE", "TE"], steps=0)
    # extract thermo values
    KE, temp, PE = summary_(params)

    # calculate total energy
    TE = KE + PE

    # time of simulation
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]

    # plot only n steps
    if steps!=0 && steps<N
        KE = KE[1:steps]
        TE = TE[1:steps]
        PE = PE[1:steps]
        t = t[1:steps]
    end

    # init figure
    fig = plot()

    # plot energies
    if "KE" in E
        plot!(t, [ustrip(i) for i in KE], label="Kinetic Energy")
    end
    if "PE" in E
        plot!(t, [ustrip(i) for i in PE], label="Potential Energy")
    end
    if "TE" in E
        plot!(t, [ustrip(i) for i in TE], label="Total Energy")
    end

    # labels
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Energy ($(UNITS.energy))")

    # return figure
    fig
end

function plot_temperature(res::MDBase.ODESolution, params; steps=0)

    # extract thermo values
    KE, temp, PE = summary_(params)

    # calculate time of simulation
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]

    # plot only n steps
    if steps!=0 && steps<N
        temp = temp[1:steps]
        t = t[1:steps]
    end

    # plot
    fig = plot(t, [ustrip(i) for i in temp], label="Temperature")

    # label
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Temperature ($(UNITS.temperature))")

    # return figure
    fig
end

function plot_temperature(res::MDBase.ODESolution, params, window)

    # set rolling window
    rw = convert(Int64, floor(window/2))

    # extract thermo values
    KE, temp, PE = summary_(params)

    # calculate time of simulation
    N = length(KE)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]

    # plot
    fig = plot(t, [ustrip(mean(temp[max(i-rw,1):min(i+rw,N)])) for i in 1:N], label="Temperature")

    # label
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Temperature ($(UNITS.temperature))")

    # return figure
    fig
end


function plot_momentum(res::MDBase.ODESolution, params)

    # get momentum (x, y, z)
    p = get_momentum_all(res, params)

    #calculate time of simulation
    N = length(res.t)
    t = res.t[1]:res.t[end]/(N-1):res.t[end]

    #plot
    fig = plot(t, [ustrip(i) for i in p[1,:]], label="P₁")
    plot!(t, [ustrip(i) for i in p[2,:]], label="P₂")
    plot!(t, [ustrip(i) for i in p[3,:]], label="P₃")

    #label
    xlabel!( "Time ($(UNITS.time))")
    ylabel!( "Momentum ($(UNITS.velocity*UNITS.mass))")

    # return figure
    fig
end
