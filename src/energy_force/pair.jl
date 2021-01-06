export acceleration!, potential_energy

# acceleration function for PairPotentials
function MDBase.acceleration!(dv, v, u, pots::Array{T,1}, params) where T <: Union{PairPotential, MLPairPotential}
    neighs = params.S.others.neighs
    BC = params.S.sim.boundary_condition
    a_ids = params.S.sim.a_ids
    @inbounds for i in 1:size(u, 2)
        ax, ay, az = u[1,i], u[2,i], u[3,i]
        N = neighs[1,i]+1
        for j in 2:N
            k = neighs[j,i]
            bx, by, bz = u[1,k], u[2,k], u[3,k]
            r, r2, dr = distance((ax, ay,az), (bx, by, bz), BC)
            mag = 0.0
            for pot in pots
                mag += potential_force(r, pot, (a_ids[i],a_ids[k]))
            end
            dv[1,i] += dr[1]*mag/(params.S.sim.mass[i]*r)
            dv[2,i] += dr[2]*mag/(params.S.sim.mass[i]*r)
            dv[3,i] += dr[3]*mag/(params.S.sim.mass[i]*r)
        end
    end
end


# energy function for PairPotentials
function MDBase.potential_energy(v, u, pots::Array{T,1}, params) where T <: Union{PairPotential, MLPairPotential}
    e = 0.0
    neighs = params.S.others.neighs
    a_ids = params.S.sim.a_ids
    @inbounds for i in 1:size(u, 2)
        ax, ay, az = u[1,i], u[2,i], u[3,i]
        N = neighs[1,i]+1
        for j in 2:N
            k = neighs[j,i]
            bx, by, bz = u[1,k], u[2,k], u[3,k]
            r, r2, dr = distance((ax, ay,az), (bx, by, bz), params.S.sim.boundary_condition)
            mag = 0.0
            for pot in pots
                mag += potential_energy(r, pot, (a_ids[i], a_ids[k]))
            end
            e += mag
        end
    end
    return e/2
end

function Base.show(stream::IO, pp::T) where T<: PairPotential
    type_ = split(split(string(T), "{")[1],".")[end]
    println(stream, type_, ":")
    println(stream, "\t",potential_expression(pp), " where")
    function preetify(symb)
        obj = getproperty(pp, symb)
        print(stream, "\t",symb,"(",unit(obj),"):\t")
        show(stream, 1ustrip(obj))
        println(stream)
    end
    for i in fieldnames(typeof(pp))
        preetify(i)
    end
end

function Docs.getdoc(t::T) where T <: MDSimulator.PairPotential
    type_ = split(split(string(T), "{")[1],".")[end]
    print("PairPotential of type $(type_).
Potential energy function is given as:
    $(potential_expression(t))

Objects of type $(type_) or in general subtype of PairPotential contains variable to hold pair potential parameters.
Please have a look below to see parameters of $(type_)
")
end
