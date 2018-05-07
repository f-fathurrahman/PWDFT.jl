using PWDFT

function test_main(kT::Float64)
    
    Nstates = 8
    Nelectrons = 6.0
    Nkpt = 2

    evals = Array{Float64}(Nstates,Nkpt)
    evals[:,1] = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    evals[:,2] = [-2.4, -1.0, -0.5, -0.2, -0.18, -0.11, -0.10, -0.05]

    Focc = Array{Float64}(Nstates,Nkpt)

    wk = zeros(Nkpt)
    wk[:] = 1.0/Nkpt

    println("\nkT = ", kT)
    spinpol = false
    println("\nspinpol = ", spinpol)
    
    Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT )
    @printf("E_fermi = %18.10f\n", E_fermi)

    for ik = 1:Nkpt
        @printf("\n")
        for ist = 1:Nstates
            @printf("[ik=%3d,ist=%3d]: %18.10f %18.10f\n", ik, ist, evals[ist,ik], Focc[ist,ik])
        end
    end
        
    #integFocc = sum_upto_E_fermi( Focc[:,ik], evals[:,ik], E_fermi[ik] )    
    #@printf("integFocc = %18.10f\n", integFocc)
    
    @printf("sum(Focc) = %18.10f\n", sum(Focc)/Nkpt)
    
    @printf("Entropy (-TS) = %18.10f\n", calc_entropy( Focc, wk, kT, is_spinpol=spinpol ))

end

test_main(0.01)

