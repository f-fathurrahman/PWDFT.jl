include("../src/smear_FD.jl")
include("../src/calc_entropy.jl")

function test_main(kT::Float64)
    Nstates = 5
    ev = Array{Float64}(Nstates)
    ev = [-2.4, -1.0, -0.5, -0.2, -0.19]
    efermi = 0.5*(ev[3] + ev[4])
    println("\nkT = ", kT)
    for spinpol in [true,false]
        println("\nspinpol = ", spinpol)
        Focc = smear_FD(ev, efermi, kT, is_spinpol=spinpol)
        for ist = 1:Nstates
            @printf("%18.10f %18.10f\n", ev[ist], Focc[ist])
        end
        @printf("sum(Focc) = %18.10f\n", sum(Focc))
        @printf("Entropy (-TS) = %18.10f\n", calc_entropy(Focc, kT, is_spinpol=spinpol))
    end
end


test_main(0.001)
test_main(0.01)
test_main(0.1)
