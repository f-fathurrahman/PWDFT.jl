function my_scf_potmix!( Ham::Hamiltonian; NiterMax=150, betamix=0.2, etot_conv_thr=1e-6 )

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

    @printf("\n")
    @printf("SCF iteration starts (with potential mixing), betamix = %f\n", betamix)
    @printf("\n")

    for iterSCF = 1:NiterMax
        
        evals = diag_LOBPCG!( Ham, psiks )
        
        Rhoe = calc_rhoe( Ham, psiks )

        # Save the old (input) potential
        Vxc_inp = copy(Ham.potentials.XC)
        VHa_inp = copy(Ham.potentials.Hartree)

        # Update potentials
        update!(Ham, Rhoe)

        # Now Ham.potentials contains new (output) potential
        
        # Calculate energies
        Ham.energies = calc_energies(Ham, psiks)

        Etot = sum(Ham.energies)

        diffEtot = abs(Etot - Etot_old)

        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        if diffEtot <= etot_conv_thr
            @printf("SCF is converged in %d iterations\n", iterSCF)
            return
        end
        Etot_old = Etot

        # Mix potentials (only Hartree and XC)
        Ham.potentials.Hartree = betamix*Ham.potentials.Hartree + (1-betamix)*VHa_inp
        Ham.potentials.XC = betamix*Ham.potentials.XC + (1-betamix)*Vxc_inp

        # Don't forget to update the total local potential
        for ispin = 1:Nspin
            for ip = 1:Npoints
                Ham.potentials.Total[ip,ispin] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                                 Ham.potentials.XC[ip,ispin]  
            end
        end


        flush(stdout)
    end
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end