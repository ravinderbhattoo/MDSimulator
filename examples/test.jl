f1(x, i) = i%2==0 ? x^3 + x^2 + x + 1/x : x
@inline f1_in(x, i) = i%2==0 ? x^3 + x^2 + x + 1/x : x

f2(x, i) = i%2==0 ? x^3 + x^2 + x + 1/x : x
@inline f2_in(x, i) = i%2==0 ? x^3 + x^2 + x + 1/x : x


function f3(x)
    s = 0.0
    for i in 1:100000000
        s += f1(x, i) + f2(x, i)
    end
    return s
end

function f3_in(x)
    s = 0.0
    for i in 1:100000000
        s += f1_in(x, i) + f2_in(x, i)
    end
    return s
end



@code_native f3(2.0)
@code_native f3_in(2.0)
