export thermo_print

function thermo_print(integrator, params)
    print("Time: $(Int(round(integrator.t/integrator.dt)))×Δτ\t")
    print("KE: $(params.global_vals[1][end]) \t")
    print("PE: $(params.global_vals[3][end]) \t")
    print("Temp: $(params.global_vals[2][end]) \t")
    print("\n")
end
