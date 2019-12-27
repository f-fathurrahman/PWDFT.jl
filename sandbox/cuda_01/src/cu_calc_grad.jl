import PWDFT: calc_grad

function calc_grad( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )
    
    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt
    #
    Focc = Ham.electrons.Focc
    #
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]
    #
    grad = CuArrays.zeros( ComplexF64, Ngw, Nstates )

    H_psi = op_H( Ham, psi )
    for ist = 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
        grad[:,ist] = Focc[ist,ikspin]*grad[:,ist]
    end

    return grad

end


