# export
export SymMat, DumpArray, nyuntam, dump2disk!

struct SymMat
    x::AbstractArray{T1,N1} where {T1,N1}
    s::Int64
    function SymMat(x)
        n = length(x)
        s = convert(Int64, sqrt((8n + 1)/4) - 1/2)
        new(x, s)
    end
end

function Base.show(stream::IO, a::SymMat)
    println(stream, "Symmetric Matrix: $(a.s)x$(a.s)")
    indi = 0
    while indi<a.s
        indj = 0
        indi += 1
        while indj<a.s
            indj += 1
            if indi==4 && a.s>7
                if indj!=4
                    print(stream, "⦙", "\t")
                else
                    print(stream, "⋱", "\t")
                end
                if indj==a.s
                    indi = a.s-3
                end
            elseif indi>=indj
                if indj==4 && a.s>7
                    print(stream, "…", "\t")
                else
                    print(stream, a[indi, indj], "\t")
                end
            else
                if indj==4 && a.s>7
                    print(stream, "…", "\t")
                else
                    print(stream, "⦿", "\t")
                end
            end
            if indj==a.s println(stream, "") end
        end
    end
end

Base.length(a::SymMat) = Base.length(a.x)
Base.size(a::SymMat) = Base.size(a.x)
Base.iterate(a::SymMat) = Base.iterate(a.x)
Base.iterate(a::SymMat, state::Int64) = Base.iterate(a.x, state)

function Base.getindex(a::SymMat, n::Int64)
    j = fld(n-1, a.s) + 1
    i = n - a.s*(j-1)
    a[i, j]
end

function Base.setindex!(a::SymMat, X, n::Int64)
    j = fld(n-1, a.s) + 1
    i = n - a.s*(j-1)
    a[i, j] = X
    nothing
end

function Base.getindex(a::SymMat, i::Int64, j::Int64)
    if i<=a.s && j<=a.s
        if j==1
            return a.x[i]
        end
        if i<j
            return a.x[convert(Int64, a.s*(i-1) - (i-2)*(i-1)/2 +  j - (i-1))]
        else
            return a.x[convert(Int64, a.s*(j-1) - (j-2)*(j-1)/2 +  i - (j-1))]
        end
    else
        @error "Trying to access array A with size [$(a.s), $(a.s)] at index [$i, $j]"
    end
end

function Base.getindex(a::SymMat, i_s::UnitRange{Int64})
    Base.getindex(a, collect(i_s))
end

function Base.getindex(a::SymMat, i_s::StepRange{Int64,Int64})
    Base.getindex(a, collect(i_s))
end

function Base.getindex(a::SymMat, i_s::UnitRange{Int64}, j_s::UnitRange{Int64})
    Base.getindex(a, collect(i_s), collect(j_s))
end

function Base.getindex(a::SymMat, i_s::Array{Int64,1})
    SymMat([a[i,j] for j in i_s for i in i_s if i<=j])
end

function Base.getindex(a::SymMat, i_s::Array{Int64,1}, j_s::Array{Int64,1})
    temp = zeros(typeof(a[1,1]), length(i_s),length(j_s))
    for j in 1:length(j_s)
        for i in 1:length(i_s)
            temp[i,j] = a[i_s[i],j_s[j]]
        end
    end
    temp
end

function Base.setindex!(a::SymMat, X, i::Int64, j::Int64)
    if i<=a.s && j<=a.s
        if j==1
            a.x[i] = X
        end
        if i<j
            a.x[convert(Int64, a.s*(i-1) - (i-2)*(i-1)/2 +  j - (i-1))] = X
        else
            a.x[convert(Int64, a.s*(j-1) - (j-2)*(j-1)/2 +  i - (j-1))] = X
        end
    else
        @error "Trying to set array A with size [$(a.s), $(a.s)] at index [$i, $j]"
    end
end

macro definemath(T, x)
    fs = [:-, :+, :exp, :sin, :cos, :tan, :log, :log10, :log2, :abs, :abs2]
    for f in fs
        eval(quote
            function Base.$f(a::$T)
                $T(Base.$f.(getproperty(a, $x)))
            end
        end)
    end

    arth = [:+, :-, :*, :/, :<, :>, :<=, :>=, :^]
    for r in arth
        eval(quote
            function Base.$r(a::$T, n::Number)
                $T(Base.$r.(getproperty(a, $x), n))
            end
            function Base.$r(n::Number, a::$T)
                $T(Base.$r.(n, getproperty(a, $x)))
            end
            function Base.$r(a::$T, b::$T)
                $T(Base.$r.(getproperty(a, $x), getproperty(b, $x)))
            end
        end
        )
    end
end

@definemath SymMat :x

Unitful.unit(a::SymMat) = Unitful.unit(a.x)
Unitful.ustrip(a::SymMat) = SymMat(1Unitful.ustrip(a.x))



# =======================
# DumpArray dumps data to file and offset array

mutable struct DumpArray{T,N} <: AbstractArray{T,N}
    x::AbstractArray{T,N}
    offset::Int64
    dumpfile::String
    dumpable::Bool
end

function Base.show(stream::IO, x::DumpArray)
    s = size(x)
    print(join(s[1:end-1], "x"));print("x($(s[end]-x.offset+1):$(s[end]))"); println(" DumpArray{$(eltype(x)),$(length(s))}:")
    println(stream, "Dump file: ",x.dumpfile)
    println(stream, "Dumpable: ",x.dumpable)
    println(stream, "Offset: ",x.offset)
    display(x.x)
end

function Base.display(x::DumpArray)
    Base.show(x)
end

function DumpArray(args...; len=1000, offset=0, dumpfile="temp_dump_array", dumpable=true)
    x = zeros(args..., len)
    DumpArray(x, offset, dumpfile, dumpable)
end

function DumpArray(dtype::DataType, args...; len=1000, offset=0, dumpfile="temp_dump_array", dumpable=true)
    x = zeros(dtype, args..., len)
    DumpArray(x, offset, dumpfile, dumpable)
end

Base.length(a::DumpArray) = prod(Base.size(a))
Base.size(a::DumpArray) = (Base.size(a.x)[1:end-1]..., Base.size(a.x)[end] + a.offset)
Base.firstindex(a::DumpArray, i::Int64) = nyuntam(a)[i]

function Base.to_indices(a::DumpArray, I::Tuple)
    I_ = []
    for (i,ind) in zip(I, 1:length(I))
        if typeof(i)==Colon
            push!(I_, nyuntam(a)[ind]:lastindex(a, ind))
        else
            push!(I_, i)
        end
    end
    return (I_...,)
end

function Base.axes(A::DumpArray)
    Base.@_inline_meta
    (map(Base.OneTo, Base.size(A)[1:end-1])..., nyuntam(A)[end]:Base.size(A)[end])
end



nyuntam(a::DumpArray) = ((1 for i in 1:length(size(a))-1)..., a.offset+1)
unavailabel(a::DumpArray) = (size(a)[1:end-1]..., a.offset)
Base.iterate(a::DumpArray) = Base.iterate(a.x)
Base.iterate(a::DumpArray, state::Int64) = Base.iterate(a.x, state)

function Base.getindex(a::DumpArray, n::Int64)
    un = unavailabel(a)
    p = prod(un)
    if n>p
        n = n - p
        Base.getindex(a.x, n)
    else
        error("Index not available. It has been dumped( <= $(unavailabel(a)) or <= $p).")
    end
end

function Base.setindex!(a::DumpArray, X, n::Int64)
    un = unavailabel(a)
    p = prod(un)
    if n>p
        n = n - p
        Base.setindex!(a.x, X, n)
    else
        error("Index not available. It has been dumped( <= $(unavailabel(a)) or <= $p).")
    end
end


function Base.getindex(a::DumpArray, I::Int64...)
    s = size(a)
    nyun = nyuntam(a)
    bool = true
    for (i, j, k) in zip(I, s, nyun)
        bool = bool && i>=k && i<=j
    end
    if bool
        return Base.getindex(a.x, I[1:end-1]..., I[end]-a.offset)
    else

        @error "Trying to access array A with size [$(prod("$i:$j, " for (i,j) in zip(nyun, s)))] at index [$(prod("$i, " for i in I))]"
    end
end

function Base.setindex!(a::DumpArray, X, I::Int64...)
    s = size(a)[end]
    nyun = nyuntam(a)[end]
    if nyun<=I[end]<=s
        return Base.setindex!(a.x, X, I[1:end-1]..., I[end] - a.offset)
    else
        @error "Trying to access array A with [$(s)] on last axis at index [$(I[end])]"
    end
end

function Base.checkbounds(a::DumpArray{<:Any,N}, I::Union{UnitRange,Int64}...) where N
    if typeof(I[end])==Int64 && I[end]==size(a)[end]+1
        true
    end
end


# function dump2disk!(a::DumpArray)
#     if a.dumpable
#         nyun = nyuntam(a)
#         s = size(a)
#         filename = "$(a.dumpfile)_($(nyun[1]):$(s[1]),$(nyun[2]):$(s[2]))"
#         @warn "Saving to file is not implemented yet."
#         a.x[:,:] .= 0
#         a.offset += size(a.x)[a.dims]
#         nothing
#     else
#         error("Cannot dump data. Please set $(typeof(a)).dumpable true.")
#     end
# end
#
#
#
function typeheirarchy(a, n::Int64)
    dtype = typeof(a)
    typeheirarchy(dtype, n)
end

function typeheirarchy(dtype::DataType, n::Int64)
    print(dtype)
    for i in 1:n
        print(" <: ")
        dtype = supertype(dtype)
        print(dtype)
        if dtype==Any
            break
        end
    end
end



#
