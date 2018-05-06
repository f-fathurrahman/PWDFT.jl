include("../src/smear_FD.jl")
include("../src/calc_Focc.jl")
include("../src/calc_entropy.jl")
include("../src/sum_upto_E_fermi.jl")

function test_main(kT::Float64)
    Nstates = 8
    Nelectrons = 6.0
    Nkpt = 1
    evals = Array{Float64}(Nstates)
    evals = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    println("\nkT = ", kT)
    spinpol = false
    println("\nspinpol = ", spinpol)
    Focc, E_fermi = calc_Focc(evals, Nelectrons, kT)
    for ist = 1:Nstates
        @printf("%18.10f %18.10f\n", evals[ist], Focc[ist])
    end
    @printf("E_fermi = %18.10f\n", E_fermi)
    integFocc = sum_upto_E_fermi( Focc, evals, E_fermi )
    @printf("integFocc = %18.10f\n", integFocc)
    @printf("sum(Focc) = %18.10f\n", sum(Focc))
    @printf("Entropy (-TS) = %18.10f\n", calc_entropy(Focc, kT, is_spinpol=spinpol))
end

test_main(0.01)
