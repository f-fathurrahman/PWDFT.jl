function my_scf!( Ham::Hamiltonian; NiterMax=150, betamix=0.2, etot_conv_thr=1e-6 )

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nspin*Nkpt
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    pspots = Ham.pspots

    psiks = rand_BlochWavefunc(Ham)
    
    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_new = zeros(Float64,Npoints,Nspin)

    Rhoe = calc_rhoe(Ham,psiks)
    update!(Ham, Rhoe)

    Ham.energies.NN = calc_E_NN(atoms)
    Ham.energies.PspCore = calc_PspCore_ene(atoms, pspots)

    evals = zeros(Nstates,Nkspin)

    Etot_old = 0.0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    for iterSCF = 1:NiterMax
        evals = diag_LOBPCG!( Ham, psiks )
        Rhoe_new = calc_rhoe( Ham, psiks )
        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        update!(Ham, Rhoe)
        Ham.energies = calc_energies(Ham, psiks)
        Etot = sum(Ham.energies)
        diffEtot = abs(Etot - Etot_old)
        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        if diffEtot <= etot_conv_thr
            @printf("SCF is converged in %d iterations\n", iterSCF)
            return
        end
        Etot_old = Etot
        flush(stdout)
    end
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end