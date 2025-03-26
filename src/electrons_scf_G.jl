# Debug version
# Similar to electrons_scf, but mixing is done in G-space
function electrons_scf_G!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.5,
    etot_conv_thr=5e-7,
    ethr_evals_last=1e-13,
    use_smearing=false,
    kT::Float64=1e-3,
    startingrhoe::Symbol=:gaussian,
    restart::Bool=false,
    print_final_ebands::Bool=true,
    starting_magnetization=nothing
)

# FIXME: use_smearing is read from electrons.use_smearing ?

    # Prepare for SCF
    # Also calculate some energy terms
    # We don't use Ham.energies to save these terms
    if (startingrhoe == :gaussian) && !restart
        Ehartree, Exc = _prepare_scf!(Ham, psiks, starting_magnetization=starting_magnetization)
    elseif (startingrhoe == :none) || restart
        Ehartree = Ham.energies.Hartree
        Exc = Ham.energies.XC
    end

    ok_paw = any(Ham.pspotNL.are_paw)

    # Reference EHxc_paw here
    if ok_paw
        EHxc_paw = Ham.pspotNL.paw.EHxc_paw
    end

    # Ham.energies.NN is only used to save this term
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    Etot_old = 0.0
    Nconv = 0

    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    dVol = CellVolume/Npoints
    Nstates = Ham.electrons.Nstates # not used?
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Vhartree = Ham.potentials.Hartree
    Vxc = Ham.potentials.XC
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk
    Nelectrons = Ham.electrons.Nelectrons
    evals = Ham.electrons.ebands

    Rhoe = Ham.rhoe
    if ok_paw
        becsum = Ham.pspotNL.paw.becsum
    end

    RhoeG = _rhoeG_from_rhoe(Ham, Rhoe) # get Rhoe in G-space
    RhoeG_in = zeros(ComplexF64, Npoints, Nspin)

    diffRhoe = 0.0

    # Mix directly in R-space
    #mixer = BroydenMixer(Rhoe, betamix, mixdim=8)
    #mixer = AdaptiveLinearMixer(Rhoe, 0.1, betamax=0.5) # not working


    Vin = zeros(Float64, Npoints, Nspin)
    Rhoe_in = zeros(Float64, Npoints, Nspin)

    if ok_paw
        becsum_in = zeros(Float64, size(Ham.pspotNL.paw.becsum))
        mixer = BroydenMixer_G(RhoeG, becsum, betamix, mixdim=8)
    else
        # Try mix in G-space
        mixer = BroydenMixer_G(RhoeG, betamix, mixdim=8)
    end

    ethr = 1e-5 # default

    is_converged = false

    # Other energy terms
    Eband = 0.0
    deband = 0.0
    descf = 0.0
    Etot = 0.0
    mTS = 0.0

    xc_calc = Ham.xc_calc
    if xc_calc.family == :metaGGA
        KEdens = zeros(Float64, Npoints, Nspin)
        KEdens_in = zeros(Float64, Npoints, Nspin)  # ????
    end

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    # additional info for spinpol calculation
    if Nspin == 2
        @printf("Initial integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Initial integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Initial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
        println("")
    end

    # Print header
    @printf("----------------------------------------------------------\n")
    @printf("     iterSCF           E              ΔE           Δρ\n")
    @printf("----------------------------------------------------------\n")

    for iterSCF in 1:NiterMax

        #println("")
        #println("Begin iterSCF = ", iterSCF)
        #println("--------------------------")
        
        # Save input/old potential
        Vin .= Vhartree .+ Vxc
        deband_hwf = -sum(Vin .* Rhoe)*dVol
        if ok_paw
            deband_hwf -= sum(Ham.pspotNL.paw.ddd_paw .* becsum)
        end
        #
        if xc_calc.family == :metaGGA
            # this is not efficient as it recalculates
            calc_KEdens!(Ham, psiks, KEdens)
            deband_hwf -= sum(xc_calc.Vtau .* KEdens[:,1])*dVol
        end

        #
        # Copy input Rhoe and RhoeG
        #
        @views Rhoe_in[:,:] .= Rhoe[:,:] # need views here?
        @views RhoeG_in[:,:] .= RhoeG[:,:]
        # also becsum if using PAW
        if ok_paw
            becsum_in[:,:,:] .= Ham.pspotNL.paw.becsum[:,:,:]
        end
        # XXX: also KEdens_in in case of metaGGA

        #
        # Diagonalization step
        #
        ethr = _calc_diag_ethr_conv(iterSCF, ethr, ethr_evals_last)
        #
        #println("\niterSCF = ", iterSCF)
        #println("Davidson diagonalization with ethr = ", ethr)
        evals[:,:] .= diag_davidson_qe!( Ham, psiks, tol=ethr )

        # Update occupation numbers from computed eigenvalues
        if use_smearing
            Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            mTS = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.mTS = mTS
            Ham.electrons.Focc = copy(Focc)
        end

        #
        # Calculate electron density and band energy
        #
        Eband = _calc_Eband(wk, Focc, evals)
        #
        calc_rhoe!( Ham, psiks, Rhoe )
        _rhoeG_from_rhoe!(Ham, Rhoe, RhoeG) # also RhoeG new
        # In case of PAW becsum is also calculated/updated here

        #println("integ Rhoe before mixing (after calc_rhoe)= ", sum(Rhoe)*dVol)
        #if Nspin == 2
        #    println("Integ magn before mixing = ", sum(Rhoe[:,1]-Rhoe[:,2])*dVol)
        #end

        # This is not used later?
        hwf_energy = Eband + deband_hwf + Ehartree + Exc + Ham.energies.NN + mTS
        if ok_paw
            hwf_energy += Ham.pspotNL.paw.EHxc_paw
        end

        # Calculate deband (using new Rhoe)
        Vin .= Vhartree .+ Vxc
        deband = -sum(Vin .* Rhoe)*dVol # TODO: use dot instead?
        if ok_paw
            deband -= sum(Ham.pspotNL.paw.ddd_paw .* becsum)
        end

        #
        if xc_calc.family == :metaGGA
            #XXX this is not efficient as it recalculates
            calc_KEdens!(Ham, psiks, KEdens)
            @assert Nspin == 1
            @views deband -= sum(xc_calc.Vtau .* KEdens[:,1])*dVol
        end

        #
        # Mix the density
        #
        #diffRhoe = norm(Rhoe - Rhoe_in)
        diffRhoe = norm(Rhoe - Rhoe_in)
        #@info "norm diffRhoe before mix = $(diffRhoe)"

        if ok_paw
            #@info "becsum before mix = $(sum(becsum))"
            do_mix!(mixer, Ham.pw, RhoeG, RhoeG_in, iterSCF,
                bec_in=becsum, bec_out=becsum_in)
            # is becsum updated here?
            #@info "becsum after mix = $(sum(becsum))"
        else
            do_mix!(mixer, Ham.pw, RhoeG, RhoeG_in, iterSCF)
        end
        _rhoe_from_rhoeG!(Ham, RhoeG, Rhoe)

        #@info "integ Rhoe after mixing = $(sum(Rhoe)*dVol)"

        #diffRhoe = dot(Rhoe - Rhoe_in, Rhoe - Rhoe_in)
        #@info "diffRhoe after mix = $(diffRhoe)"

        # Check convergence here? (using diffRhoe)

        # XXX Rhoe and becsum are used here
        Ehartree, Exc = update_from_rhoe!(Ham, psiks, Rhoe)

        descf = -sum( (Rhoe_in .- Rhoe).*(Vhartree .+ Vxc) )*dVol
        # XXX: metagga contribution is not included in descf yet !!!
        if ok_paw
            descf -= sum(Ham.pspotNL.paw.ddd_paw .* (becsum_in - becsum))
            #@info "sum diff becsum $(sum(becsum_in - becsum))"
        end

        Etot = Eband + deband + Ehartree + Exc + Ham.energies.NN + descf + mTS
        if ok_paw
            Etot += Ham.pspotNL.paw.EHxc_paw
        end
    
        #
        diffEtot = abs(Etot - Etot_old)
        @printf("SCF_G: %5d  %18.10f  %12.5e  %12.5e\n", iterSCF, Etot, diffEtot, diffRhoe)
        if Nspin == 2
            println("integ magn = ", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
        end
        if mixer.is_converged
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2 # FIXME: make this a parameter
            @printf("SCF_G: Total energy is converged in %d iterations\n", iterSCF)
            is_converged = true
            break
        end

        Etot_old = Etot
        flush(stdout)
    end

    println()
    println(">>>> Final result:")
    println()
    println("-----------------------")
    println("Energy components in Ry")
    println("-----------------------")
    @printf("Eband    = %18.10f Ry\n", Eband*2)
    @printf("deband   = %18.10f Ry\n", deband*2)
    @printf("descf    = %18.10f Ry\n", descf*2)
    @printf("-----------------------------\n")
    @printf("OneEle   = %18.10f Ry\n", 2*(Eband + deband))
    @printf("Ehartree = %18.10f Ry\n", 2*Ehartree)
    @printf("Exc      = %18.10f Ry\n", 2*Exc)
    @printf("NN       = %18.10f Ry\n", 2*Ham.energies.NN)
    @printf("mTS      = %18.10f Ry\n", 2*mTS)
    if ok_paw
        @printf("EHxc_paw = %18.10f Ry\n", 2*Ham.pspotNL.paw.EHxc_paw)
    end
    @printf("-----------------------------\n")
    @printf("! Total  = %18.10f Ry\n", 2*Etot)

    if ok_paw
        @printf("Total all electrons = %18.10f Ry\n", 2*(Etot + Ham.pspotNL.paw.total_core_energy))
    end

    if Nspin == 2
        println("integ magn = ", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
    end

    # TODO
    # Also print the Kohn-Sham orbital energies using similar format
    # as in pw.x


    if !is_converged
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    # Compare the energy using the usual formula (not using double-counting)
    calc_energies!(Ham, psiks)
    if use_smearing
        Ham.energies.mTS = mTS
    end
    println("\nUsing original formula for total energy")
    println(Ham.energies, use_smearing=use_smearing, is_paw=ok_paw)
    
    if print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")
    end

    return
end


