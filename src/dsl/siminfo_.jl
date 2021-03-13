macro litstr(name, func)
    name_ = name*"_str"
    ex = Meta.parse(
        """macro $name_(s)
            lookup = create_lookup(s)
            dict = Dict{Symbol, Number}()
            for (k, v) in lookup
                if string(k)[1]!="_"
                    dict[k] = v
                end
            end
            return :($func(;\$(dict)...))
        end""")
    eval(ex)
end

function resolve_expr!(ex, lookup)
    for i in 2:length(ex.args)
        if typeof(ex.args[i])==Symbol
            ex.args[i] = lookup[ex.args[i]]
        elseif typeof(ex.args[i])==Expr
            resolve_expr!(ex.args[i], lookup)
        end
    end
end

function create_lookup(s)
    lookup = Dict{Symbol, Union{Number,Array{<:Number},Symbol,Expr}}()
    if s!=""
        params = split(s, ";")
        for p in params
            ex = Meta.parse(p)
            lookup[ex.args[1]] = ex.args[2]
        end
    end

    for (k, v) in lookup
        if typeof(v)==Symbol
            lookup[k] = lookup[v]
        end
        if typeof(v)==Expr
            resolve_expr!(v, lookup)
            lookup[k] = eval(lookup[k])
        end
    end
    return lookup
end
