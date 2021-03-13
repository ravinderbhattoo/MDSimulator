function grd(f, x)
    y, back = Zygote.pullback(f, x)
    back(1)
end

function jac(f, x)
    y = f(x)
    etype_ = eltype(promote(y[1], x[1])[1])
    # G = Zygote.Buffer([etype_(1)], length(y), length(x))
    A = zeros(length(y), length(x))
    for i in 1:length(y)
        @show i
        y_, back = Zygote.pullback((x) -> f(x)[i], x)
        @show back(1)
        A[i,:] .= reshape(back(1)[1], (:))
    end
    A
end
