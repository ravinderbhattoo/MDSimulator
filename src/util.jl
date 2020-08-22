function summary(res)
    KE, temp, PE = res.prob.p[end][end]
    try
        return KE, temp, PE, KE+PE
    finally
        return KE, temp, PE
    end
end
