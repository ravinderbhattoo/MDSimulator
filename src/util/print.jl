export thermo_print, thermo_headers, set_thermo, describe

function set_thermo()
    set_thermo(Symbol[:temp, :ke, :pe])
end

function set_thermo_(keys_::Array{Symbol,1}; ckeys_::Array{String,1}=String[], output="none", sigdigits=6)
    key_loc = Dict{String, String}()
    for i in keys_
        key_loc[string(i)] = "a.$(i)[n]"
    end
    for i in ckeys_
        key_loc[i] = "a.custom[\"$(i)\"][n]"
    end

    ks = collect(keys(key_loc))
    sort!(ks)

    if output=="none"
        p3 = "print(\"n=\", n, '\t', '\t'); " * mapreduce((x)->"print(\"$x\", '=', round($(key_loc[x]), sigdigits=$sigdigits), '\t','\t'); ", *,ks)
    elseif output=="csv"
        p3 = "print(n, ','); " * mapreduce((x)->"print(round($(key_loc[x]), sigdigits=$sigdigits), ','); ", *,ks)
    elseif output=="dtab-delim"
        p3 = "print(n, '\t', '\t'); " * mapreduce((x)->"print(round($(key_loc[x]), sigdigits=$sigdigits), '\t','\t'); ", *,ks)
    end

    p1 = "(a, n) -> begin "
    p2 = ""
    p4 = "println(); end"
    return eval(Meta.parse(p1*p2*p3*p4))
end


function set_thermo_headers_(keys_::Array{Symbol,1}, ckeys_::Array{String,1}; output="none")
    key_ = String[]
    for i in keys_
        push!(key_, string(i))
    end
    for i in ckeys_
        push!(key_, string(i))
    end
    sort!(key_)

    if output=="none"
        p3 = ""
    elseif output=="csv"
        p3 = "print(\"n,\"); " * mapreduce((x)->"print(\"$x,\"); ", *,key_)
    elseif output=="dtab-delim"
        p3 = "print(\"n\t\t\"); " * mapreduce((x)->"print(\"$x\t\t\"); ", *,key_)
    end

    p1 = "() -> begin "
    p2 = "println(\"\nThermo output:\"); "
    p4 = "println(); end"
    return eval(Meta.parse(p1*p2*p3*p4))
end

function thermo_print end
function thermo_headers end

function set_thermo(keys_::Array{Symbol,1}; ckeys_::Array{String,1}=String[], kwargs...)
    local f = set_thermo_(keys_; ckeys_=ckeys_, kwargs...)
    local g = set_thermo_headers_(keys_, ckeys_; kwargs...)
    eval(:(MDSimulator.thermo_print(x, n) = $f(x, n)))
    eval(:(MDSimulator.thermo_headers() = $g()))
end

set_thermo([:temp])


function describe(sol::T, params) where T<:MDBase.ODESolution
    delta = params.S.sim.save_every
    n = length(sol.t)-1
    N = length(params.S.others.thermo_vals.ke)-1
    if sol.retcode == :Unstable
        println("Simulation is unstable.")
    end
    if sol.retcode == :Success
        println("Simulation is successful.")
    end
    println("Number of runs: $(n*delta)/$(N*delta)")
end
