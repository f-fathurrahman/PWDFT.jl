function calc_grad( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    
    ik = Ham.ik
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc
    pw = Ham.pw
    #
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]
    Ω = pw.Ω
    Ns = pw.Ns
    #
    grad = zeros( Complex128, Ngw, Nstates )

    H_psi = op_H( Ham, psi )
    for ist = 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
        grad[:,ist] = Focc[ist,ik]*grad[:,ist]
    end
    return grad

end
