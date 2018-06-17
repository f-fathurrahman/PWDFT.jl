 function calc_grad( Ham::PWHamiltonian, psi::Array{ComplexF64,2} )
    
    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt
    #
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc
    pw = Ham.pw
    #
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]
    Ω = pw.Ω
    Ns = pw.Ns
    #
    grad = zeros( ComplexF64, Ngw, Nstates )

    H_psi = op_H( Ham, psi )
    for ist = 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
        grad[:,ist] = Focc[ist,ikspin]*grad[:,ist]
    end

    # the usual case of constant occupation numbers
    if all(Focc[:,ikspin] .== 2.0) == true || all(Focc[:,ikspin] .== 1.0) == true
        # immediate return
        return grad
    end

    # Calculate reduced Hamiltonian
    Hred = zeros(ComplexF64,Nstates,Nstates)
    for ist = 1:Nstates
        for jst = ist:Nstates
            Hred[ist,jst] = dot(psi[:,jst], H_psi[:,ist])
            Hred[jst,ist] = Hred[ist,jst]
        end
    end

    # Additional contribution
    for ist = 1:Nstates
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] + 0.5*Hred[ist,jst]*
                          (Focc[ist,ikspin]-Focc[jst,ikspin])*psi[:,jst]
        end
    end

    return grad

end
