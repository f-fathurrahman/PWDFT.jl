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

    println("\nkT = ", kT)
    
    Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT )
    @printf("E_fermi = %18.10f\n", E_fermi)

    for ik = 1:Nkpt
        @printf("\n")
        for ist = 1:Nstates
            @printf("[ik=%3d,ist=%3d]: %18.10f %18.10f\n", ik, ist, evals[ist,ik], Focc[ist,ik])
        end
    end
        
    @printf("sum(Focc) = %18.10f\n", sum(Focc)/Nkpt)
    @printf("Entropy (-TS) = %18.10e\n", calc_entropy( Focc, wk, kT, Nspin=Nspin ))

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

test_no_spin(0.01)
test_spin(0.01)

