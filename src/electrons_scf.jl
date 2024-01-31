# FIXME: To be included in the main package file
include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("smooth_to_dense.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")
include("ortho_with_S.jl")
include("diag_davidson_qe_v2.jl")

# An implementation of electrons_scf subroutine in PWSCF
# Total energy is calculated using double-counting formula
function electrons_scf!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.5,
    etot_conv_thr=1e-6,
    ethr_evals_last=1e-13,
    use_smearing=false,
    kT::Float64=1e-3,
    startingrhoe::Symbol=:gaussian,
    restart::Bool=false
)

    # Prepare for SCF
    # Also calculate some energy terms
    # We don't use Ham.energies to save these terms
    if (startingrhoe == :gaussian) && !restart
        Ehartree, Exc = _prepare_scf!(Ham, psiks)
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
    Nstates = Ham.electrons.Nstates

    Rhoe = Ham.rhoe
    println("Initial integ Rhoe = ", sum(Rhoe)*dVol)

    diffRhoe = 0.0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    Vin = zeros(Float64, Npoints)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Vhartree = Ham.potentials.Hartree
    Vxc = Ham.potentials.XC
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk
    Nelectrons = Ham.electrons.Nelectrons

    evals = Ham.electrons.ebands

    # Mix directly in R-space
    mixer = BroydenMixer(Rhoe, betamix, mixdim=8)
    #mixer = AdaptiveLinearMixer(Rhoe, 0.1, betamax=0.5) # not working

    Rhoe_in = zeros(Float64, size(Rhoe))
    if ok_paw
        becsum_in = zeros(Float64, size(Ham.pspotNL.becsum))
        mixer_becsum = BroydenMixer(Ham.pspotNL.becsum, betamix, mixdim=8)
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


    # Print header
    @printf("----------------------------------------------------------\n")
    @printf("     iterSCF           E              ΔE           Δρ\n")
    @printf("----------------------------------------------------------\n")

    for iterSCF in 1:NiterMax
        
        @views Vin[:] .= Vhartree[:] + Vxc[:,1]
        @views deband_hwf = -sum(Vin .* Rhoe[:,1])*dVol
        if ok_paw
            ddd_paw = Ham.pspotNL.paw.ddd_paw
            becsum = Ham.pspotNL.becsum
            deband_hwf -= sum(ddd_paw .* becsum)
        end
        #
        if xc_calc.family == :metaGGA
            # this is not efficient as it recalculates
            calc_KEdens!(Ham, psiks, KEdens)
            deband_hwf -= sum(xc_calc.Vtau .* KEdens[:,1])*dVol
        end

        # Copy input rhoe
        @views Rhoe_in[:,:] .= Rhoe[:,:]
        if ok_paw
            becsum_in[:,:,:] .= Ham.pspotNL.becsum[:,:,:]
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

        if use_smearing
            Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            mTS = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        #
        # Calculate electron density and band energy
        #
        Eband = _calc_Eband(wk, Focc, evals)
        #Rhoe[:,:] = calc_rhoe( Ham, psiks )
        calc_rhoe!( Ham, psiks, Rhoe )
        # In case of PAW becsum is also calculated/updated here

        # This is not used later?
        hwf_energy = Eband + deband_hwf + Ehartree + Exc + Ham.energies.NN + mTS
        if ok_paw
            hwf_energy += Ham.pspotNL.paw.EHxc_paw
        end

        # Calculate deband (using new Rhoe)
        @views Vin[:] .= Vhartree[:] + Vxc[:,1]
        @views deband = -sum(Vin .* Rhoe[:,1])*dVol # TODO: use dot instead?
        if ok_paw
            ddd_paw = Ham.pspotNL.paw.ddd_paw
            becsum = Ham.pspotNL.becsum
            deband -= sum(ddd_paw .* becsum)
        end

        #
        if xc_calc.family == :metaGGA
            # this is not efficient as it recalculates
            calc_KEdens!(Ham, psiks, KEdens)
            @views deband -= sum(xc_calc.Vtau .* KEdens[:,1])*dVol
        end

        #
        # Mix the density
        #
        #diffRhoe = norm(Rhoe - Rhoe_in)
        diffRhoe = dot(Rhoe - Rhoe_in, Rhoe - Rhoe_in)
        
        do_mix!(mixer, Rhoe, Rhoe_in, iterSCF)

        # Check convergence here? (using diffRhoe)

        #
        Ehartree, Exc = update_from_rhoe!(Ham, psiks, Rhoe)

        @views descf = -sum( (Rhoe_in[:,1] .- Rhoe[:,1]).*(Vhartree + Vxc[:,1]) )*dVol
        # XXX: metagga contribution is not included in descf yet !!!
        if ok_paw
            ddd_paw = Ham.pspotNL.paw.ddd_paw
            becsum = Ham.pspotNL.becsum
            descf -= sum(ddd_paw .* (becsum_in - becsum))
        end

        Etot = Eband + deband + Ehartree + Exc + Ham.energies.NN + descf + mTS
        if ok_paw
            Etot += Ham.pspotNL.paw.EHxc_paw
        end
    
        #
        diffEtot = abs(Etot - Etot_old)
        @printf("SCF: %5d  %18.10f  %12.5e  %12.5e\n", iterSCF, Etot, diffEtot, diffRhoe)
        if diffEtot <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2
            @printf("SCF: Total energy is converged in %d iterations\n", iterSCF)
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
    
    @printf("\n")
    @printf("----------------------------\n")
    @printf("Final Kohn-Sham eigenvalues:\n")
    @printf("----------------------------\n")
    @printf("\n")
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")

    return
end


# Initialize Rhoe, potentials
function _prepare_scf!(Ham, psiks)
    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham)
    # Also initialize becsum in case of PAW
    if any(Ham.pspotNL.are_paw)
        PAW_atomic_becsum!(Ham)
    end

    # Update the potentials
    Ehartree, Exc = update_from_rhoe!( Ham, psiks, Rhoe, RhoeG )

    return Ehartree, Exc
end


function _calc_Eband(wk, Focc, evals)
    Nstates = size(evals, 1)
    Nkspin = size(evals, 2)
    Nkpt = size(wk, 1)
    Nspin = Int64(Nkspin/Nkpt)
    Eband = 0.0
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        for ist in 1:Nstates
            Eband = Eband + wk[ik]*Focc[ist,ikspin]*evals[ist,ikspin]
        end
    end
    return Eband
end

# determine convergence criteria for diagonalization
function _calc_diag_ethr_conv(iterSCF, ethr_current, ethr_evals_last)
    if iterSCF == 1
        ethr = 0.1
    elseif iterSCF == 2
        ethr = 0.01
    else
        ethr = ethr_current/5.0
        ethr = max( ethr, ethr_evals_last )
    end
    return ethr
end

