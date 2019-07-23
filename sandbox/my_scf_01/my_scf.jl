function my_scf!( Ham::Hamiltonian; NiterMax=150, betamix=0.2, etot_conv_thr=1e-6, log_file="scf_history.dat" )
    
    psiks = rand_BlochWavefunc(Ham)
    
    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)
    
    Rhoe_new = similar(Rhoe)
    
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    Ham.energies.PspCore = calc_PspCore_ene(Ham.atoms, Ham.pspots)

    Etot_old = 0.0
    Nconv = 0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")
    
    f = open(log_file, "w")

    for iterSCF = 1:NiterMax
        _ = diag_LOBPCG!( Ham, psiks )
        Rhoe_new = calc_rhoe( Ham, psiks )
        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        update!(Ham, Rhoe)
        Ham.energies = calc_energies(Ham, psiks)
        Etot = sum(Ham.energies)
        diffEtot = abs(Etot - Etot_old)
        @printf(f, "%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        if diffEtot <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2
            @printf("SCF is converged in %d iterations\n", iterSCF)
            close(f)
            return
        end
        Etot_old = Etot
        flush(stdout)
    end
    close(f)
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end