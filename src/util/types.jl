# export
export SymMat

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
#
