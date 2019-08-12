function grad_obj_function_v1!(
    Ham::Hamiltonian,
    psiks_::BlochWavefunc,
    g::BlochWavefunc,
    kT::Float64=1e-3,    
    skip_ortho=false
)

    psiks = copy(psiks_)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons    
    wk = Ham.pw.gvecw.kpoints.wk

    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    Rhoe_old = copy( Ham.rhoe )
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    evals = zeros(Nstates, Nkspin)

    # Rotate the subspace to calculate eigenvalues
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        evals[:,i], evecs = eigen(Hr)
    end

    Focc_old = copy( Ham.electrons.Focc )

    Ham.electrons.Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
    Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )

    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        g[i] = calc_grad_v1( Ham, psiks[i] )
    end

    update!( Ham, Rhoe_old )
    Ham.electrons.Focc = Focc_old

    return

end

function calc_grad_v1( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    
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
    grad = zeros( ComplexF64, Ngw, Nstates )

    H_psi = op_H( Ham, psi )
    for ist = 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst = 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * psi[:,jst]
        end
        grad[:,ist] = Focc[ist,ikspin]*grad[:,ist]
    end

    #F = Matrix(Diagonal(Focc[:,ikspin]))
    #ℍ = psi' * H_psi
    #HFH = ℍ*F - F*ℍ
    #ℚ = 0.5*HFH
    #grad[:,:] = grad[:,:] + psi*ℚ

    return grad
end