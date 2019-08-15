function obj_function_v1!(
    Ham::Hamiltonian,
    psiks_::BlochWavefunc;
    kT::Float64=1e-3,
    skip_ortho=false
)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons    
    wk = Ham.pw.gvecw.kpoints.wk

    psiks = copy(psiks_)

    # orthonormalize if it is needed
    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    ebands = zeros(Nstates, Nkspin)

    # Rotate the subspace to calculate eigenvalues
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        ebands[:,i], evecs = eigen(Hr)
        psiks[i] = psiks[i]*evecs
    end

    Focc_old = copy(Ham.electrons.Focc)

    Ham.electrons.Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, ebands, Nspin )
    Entropy = calc_entropy( wk, kT, ebands, E_fermi, Nspin )

    Rhoe_old = copy( Ham.rhoe )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    
    energies = calc_energies( Ham, psiks )
    energies.mTS = Entropy

    # set Rhoe back
    update!( Ham, Rhoe_old )
    #Ham.electrons.Focc = copy(Focc_old)
    Ham.electrons.ebands = ebands

    return sum( energies )

end