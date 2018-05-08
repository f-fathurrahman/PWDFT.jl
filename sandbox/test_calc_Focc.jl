using PWDFT

function test_main(kT::Float64)
    
    Nstates = 8
    Nelectrons = 6.0
    Nkpt = 2

    evals = Array{Float64}(Nstates,Nkpt)
    evals[:,1] = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    evals[:,2] = [-2.4, -1.0, -0.5, -0.2, -0.18, -0.11, -0.10, -0.05]

    Focc = Array{Float64}(Nstates,Nkpt)

    println("\nkT = ", kT)
    spinpol = false
    println("\nspinpol = ", spinpol)
    
    E_fermi = zeros(Nkpt)

    for ik = 1:Nkpt
        Focc[:,ik], E_fermi[ik] = calc_Focc( evals[:,ik], Nelectrons, kT )
        @printf("%d E_fermi = %18.10f\n", ik, E_fermi[ik])

        for ist = 1:Nstates
            @printf("%18.10f %18.10f\n", evals[ist,ik], Focc[ist,ik])
        end
        
        integFocc = sum_upto_E_fermi( Focc[:,ik], evals[:,ik], E_fermi[ik] )    
        @printf("integFocc = %18.10f\n", integFocc)
        @printf("sum(Focc) = %18.10f\n", sum(Focc[:,ik]))
        @printf("Entropy (-TS) = %18.10f\n", calc_entropy( Focc[:,ik], kT, is_spinpol=spinpol ))
    end

end

test_main(0.01)

