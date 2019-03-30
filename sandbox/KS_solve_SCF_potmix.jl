function KS_solve_SCF_potmix!(
    Ham::Hamiltonian;
    NiterMax=150,
    betamix=0.2,
    verbose=true,
    print_final_ebands=true,
    print_final_energies=true,
    savewfc=false,
    use_smearing=false,
    ETOT_CONV_THR=1e-6
)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nspin*Nkpt
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    pspots = Ham.pspots

    psiks = rand_BlochWavefunc(Ham)
    
    Rhoe = zeros(Float64,Npoints,Nspin)

    Vxc_inp = zeros(Float64,Npoints,Nspin)
    VHa_inp = zeros(Float64,Npoints)

    Rhoe = calc_rhoe(Ham,psiks)
    update!(Ham, Rhoe)

    Ham.energies.NN = calc_E_NN(atoms)
    Ham.energies.PspCore = calc_PspCore_ene(atoms, pspots)

    evals = zeros(Nstates,Nkspin)

    Etot_old = 0.0

    CONVERGED = 0

    @printf("\n")
    @printf("SCF iteration starts (with potential mixing), betamix = %f\n", betamix)
    @printf("\n")

    for iterSCF = 1:NiterMax
        
        evals = diag_LOBPCG!( Ham, psiks )
        
        Rhoe[:,:] = calc_rhoe( Ham, psiks )

        # Save the old (input) potential
        Vxc_inp[:,:] = Ham.potentials.XC
        VHa_inp[:] = Ham.potentials.Hartree

        # Update potentials
        update!(Ham, Rhoe)

        # Now Ham.potentials contains new (output) potential
        
        # Calculate energies
        Ham.energies = calc_energies(Ham, psiks)

        Etot = sum(Ham.energies)

        diffEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)

        if diffEtot < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            @printf("SCF is converged in %d iterations\n", iterSCF)
            break
        end
        Etot_old = Etot

        # Mix potentials (only Hartree and XC)
        Ham.potentials.Hartree = betamix*Ham.potentials.Hartree + (1-betamix)*VHa_inp
        Ham.potentials.XC = betamix*Ham.potentials.XC + (1-betamix)*Vxc_inp

        flush(stdout)
    end

    if CONVERGED < 2
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    Ham.electrons.ebands = evals

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)
    end

    if verbose && print_final_energies
        @printf("\n")
        @printf("-------------------------\n")
        @printf("Final Kohn-Sham energies:\n")
        @printf("-------------------------\n")
        @printf("\n")
        println(Ham.energies, use_smearing=use_smearing)
    end

    if savewfc
        for ikspin = 1:Nkpt*Nspin
            wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","w")
            write( wfc_file, psiks[ikspin] )
            close( wfc_file )
        end
    end

    return
end