using Printf
using PWDFT

function test_no_spin(kT::Float64)
    
    Nstates = 8
    Nelectrons = 8.0
    Nkpt = 2
    Nspin = 1

    evals = zeros(Float64,Nstates,Nkpt)
    evals[:,1] = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    evals[:,2] = [-2.4, -1.0, -0.5, -0.2, -0.18, -0.11, -0.10, -0.05]

    Focc = zeros(Float64,Nstates,Nkpt)

    wk = zeros(Nkpt)
    wk[:] .= 1.0/Nkpt

    Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT )

    for ik = 1:Nkpt
        @printf("\n")
        for ist = 1:Nstates
            @printf("[ik=%3d,ist=%3d]: %18.10f %18.10f\n", ik, ist, evals[ist,ik], Focc[ist,ik])
        end
    end
    
    @printf("sum(Focc) = %18.10f\n", sum(Focc)/Nkpt)

    println("wk = ", wk)

    println("\nkT = ", kT)
    @printf("E_fermi = %18.10f\n", E_fermi)

    @printf("Entropy (-TS)    = %18.10e\n", calc_entropy( Focc, wk, kT, Nspin=Nspin ))
    
    @printf("Entropy (-TS) v2 = %18.10e\n", calc_entropy_v2( wk, kT, evals, E_fermi, Nspin ))

end

function test_spin(kT::Float64)
    
    Nstates = 8
    Nelectrons = 6.0
    Nkpt = 2
    Nspin = 2
    Nkspin = Nkpt*Nspin

    evals = zeros(Float64,Nstates,Nkspin)
    # spin up
    evals[:,1] = [-2.4, -1.0, -0.5, -0.2, -0.19, -0.10, -0.05, -0.05]
    evals[:,2] = [-2.4, -1.0, -0.5, -0.2, -0.18, -0.11, -0.10, -0.05]
    # spin dn
    evals[:,3] = [-2.5, -1.1, -0.6, -0.3, -0.20, -0.11, -0.06, -0.04]
    evals[:,4] = [-2.5, -1.1, -0.6, -0.3, -0.19, -0.12, -0.10, -0.01]

    Focc = zeros(Float64,Nstates,Nkspin)

    wk = zeros(Nkpt)
    wk[:] .= 1.0/Nkpt

    println("\nkT = ", kT)
    
    Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT, Nspin=Nspin, verbose=true)
    @printf("E_fermi = %18.10f\n", E_fermi)

    integFocc = 0.0
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        @printf("\n")
        for ist = 1:Nstates
            @printf("[ispin=%d,ik=%3d,ist=%3d]: %18.10f %18.10f\n",
                     ispin, ik, ist, evals[ist,ikspin], Focc[ist,ikspin])
        end
    end
    end

    @printf("sum(Focc) = %18.10f\n", sum(Focc)/Nkpt)
    @printf("Entropy (-TS) = %18.10f\n", calc_entropy( Focc, wk, kT, evals, E_fermi, Nspin=Nspin ))

end


function test_Al_no_spin( kT::Float64 )
    
    Nstates = 6
    Nelectrons = 3.0
    Nkpt = 4
    Nspin = 1

    evals = zeros(Float64,Nstates,Nkpt)
    evals[:,1] = [-0.0418808004, 0.8381197256, 0.8388139806,
                   0.8395961149, 0.8956036622, 0.8977654851]
    evals[:,2] = [0.0675054958, 0.3866374362, 0.7694392074,
                  0.7740564059, 0.8675063924, 0.8849688613]
    evals[:,3] = [0.1034191812, 0.5290915192, 0.5945346478,
                  0.6083846984, 0.6098449963, 0.7441310954]
    evals[:,4] = [0.2442976005, 0.3257127101, 0.3872675749,
                  0.6417726442, 0.7048505310, 1.0656432986]

    Focc = zeros(Float64,Nstates,Nkpt)

    wk = zeros(Nkpt)
    wk[:] .= 1.0/Nkpt

    Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT )

    for ik = 1:Nkpt
        @printf("\n")
        for ist = 1:Nstates
            @printf("[ik=%3d,ist=%3d]: %18.10f %18.10f\n", ik, ist, evals[ist,ik], Focc[ist,ik])
        end
    end
    
    @printf("sum(Focc) = %18.10f\n", sum(Focc)/Nkpt)

    println("wk = ", wk)

    println("\nkT = ", kT)
    @printf("E_fermi = %18.10f\n", E_fermi)

    @printf("Entropy (-TS)    = %18.10e\n", calc_entropy( Focc, wk, kT, Nspin=Nspin ))
    
    @printf("Entropy (-TS) v2 = %18.10e\n", calc_entropy_v2( wk, kT, evals, E_fermi, Nspin ))

end

#test_no_spin(0.01)
#test_spin(0.01)

test_Al_no_spin(0.001)

