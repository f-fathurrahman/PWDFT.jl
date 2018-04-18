function calc_grad( Ham::PWHamiltonian, psi::Array{Complex128,2} )

    potentials = Ham.potentials
    Focc = Ham.electrons.Focc
    pw = Ham.pw
    #
    Ngwx    = size(psi)[1]
    Nstates = size(psi)[2]
    Ω = pw.Ω
    Ns = pw.Ns
    #
    grad = zeros( Complex128, Ngwx, Nstates )

    H_psi = op_H( Ham, psi )
    for ist = 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
        grad[:,ist] = Focc[ist]*grad[:,ist]
    end
    return grad

end
