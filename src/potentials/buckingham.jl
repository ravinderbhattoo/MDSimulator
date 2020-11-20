# exports
export potential_force, potential_energy, energy_unit, force_unit

######################################################
# Buckingham Potential                               #
######################################################

export BuckinghamParameters
struct BuckinghamParameters{AType <: Number, BType <: Number, CType <: Number, lType <: Number} <: PairPotential
    A::Array{AType,2}
    B::Array{BType,2}
    C::Array{CType,2}
    R::Array{lType,2}
    equillibrium::Array{lType,2}
    well_energy::Array{AType,2}
    energy_shift::Array{AType,2}
end

function BuckinghamParameters(A::Array{T1,2}, B::Array{T2,2}, C::Array{T3,2}, R::Array{T4,2}) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
    energy_shift = A.*exp.(-B.*R) .- C./R.^6
    r = R/2
    for i in 1:100
        r = (6C./(A.*B).*exp.(-B.*r)).^(1/7)
    end
    well_energy = A.*exp.(-B.*r) .- C./r.^6
    BuckinghamParameters(A, B, C, R, r, well_energy, -energy_shift)
end

function BuckinghamParameters()
    A = 1.69*1.0e-18u"kJ"*CONSTANTS.Avogadro
    B = 1/27.3u"pm"
    C = 102*1.0e-28u"kJ"*CONSTANTS.Avogadro*(1u"pm")^6
    R = 500u"pm"
    BuckinghamParameters(A, B, C, R)
end

function BuckinghamParameters(A::T1, B::T2, C::T3, R::T4) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
    ONE = ones(1,1)
    BuckinghamParameters(A*ONE, B*ONE, C*ONE, R*ONE)
end

function BuckinghamParameters(A::Array{T1,N}, B::Array{T2,N}, C::Array{T3,N}, R::Array{T4,N}) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number, N}
    l = convert(Int64, sqrt(length(A)+1))
    function fillit(a, l)
        out = zeros(l,l)
        ind = 1
        for i in 1:l
            for j in i:l
                out[i,j] = a[ind]
                out[j,i] = a[ind]
                ind += 1
            end
        end
        return out
    end
    BuckinghamParameters(fillit(A, l), fillit(B, l), fillit(C, l), fillit(R, l))
end

function potential_energy(r::Number, pot::BuckinghamParameters, pair::Tuple{Int64,Int64})
    i,j  = pair
    if r<pot.R[i,j]
        return pot.A[i,j]*exp(-pot.B[i,j]*r) - pot.C[i,j]/r^6 + pot.energy_shift[i,j]
    else
        return 0.0
    end
end

function potential_force(r::Number, pot::BuckinghamParameters, pair::Tuple{Int64,Int64})
    i,j = pair
    if r<pot.R[i,j]
        return pot.A[i,j]*pot.B[i,j]*exp(-pot.B[i,j]*r)-6pot.C[i,j]/r^7
    else
        return 0.0
    end
end

energy_unit(pot::BuckinghamParameters) = unit(pot.well_energy[1])
force_unit(pot::BuckinghamParameters) = unit(pot.C[1]/pot.R[1]^7)

function Base.show(stream::IO, pp::BuckinghamParameters)
    println(stream, "Buckingham:")
    println(stream, "\tV(r) = A e⁻ᴮʳ - C r⁻⁶ where")
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

######################################################
# Buckingham Potential                               #
######################################################


#
