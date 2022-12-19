# FIXME: To be included in the main package file
include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("smooth_to_dense.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")
include("calc_rhoe_uspp.jl")
include("ortho_with_S.jl")
include("diag_davidson_qe_v2.jl")

# An implementation of electrons_scf subroutine in PWSCF
# Total energy is calculated using double-counting formula
function electrons_scf!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.5,
    etot_conv_thr=1e-6,
    ethr_evals_last=1e-13
)

    # Prepare for SCF
    # Also calculate some energy terms
    # We don't use Ham.energies to save these terms
    Ehartree, Exc, Evtxc = _prepare_scf!(Ham, psiks)

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

    evals = Ham.electrons.ebands
        
    # Mix directly in R-space
    mixer = BroydenMixer(Rhoe, betamix, mixdim=8)

    Rhoe_in = zeros(Float64,size(Rhoe))

    ethr = 1e-5 # default

    is_converged = false

    # Other energy terms
    Eband = 0.0
    deband = 0.0
    descf = 0.0
    Etot = 0.0

    for iterSCF in 1:NiterMax
        
        @views Vin[:] .= Vhartree[:] + Vxc[:,1]
        deband_hwf = -sum(Vin .* Rhoe[:,1])*dVol

        # Copy input rhoe
        @views Rhoe_in[:,:] .= Rhoe[:,:]

        #
        # Diagonalization step
        #
        ethr = _calc_diag_ethr_conv(iterSCF, ethr, ethr_evals_last)
        #
        println("\niterSCF = ", iterSCF)
        println("Davidson diagonalization with ethr = ", ethr)
        evals[:,:] .= diag_davidson_qe!( Ham, psiks, tol=ethr )

        #
        # Calculate electron density and band energy
        #
        Eband = _calc_Eband(wk, Focc, evals)
        Rhoe[:,:] = calc_rhoe_uspp( Ham, psiks )
        println("integ output Rhoe = ", sum(Rhoe)*dVol)

        # This is not used later?
        hwf_energy = Eband + deband_hwf + Ehartree + Exc + Ham.energies.NN
        # entropy term missing

        # Calculate deband (using new Rhoe)
        @views Vin[:] .= Vhartree[:] + Vxc[:,1]
        deband = -sum(Vin .* Rhoe[:,1])*dVol

        #
        # Mix the density
        #
        #diffRhoe = norm(Rhoe - Rhoe_in)
        diffRhoe = dot(Rhoe - Rhoe_in, Rhoe - Rhoe_in)
        @printf("Before mix: diffRhoe = %e\n", diffRhoe)
        do_mix!(mixer, Rhoe, Rhoe_in, iterSCF)
        println("integ Rhoe after mix: ", sum(Rhoe)*dVol)

        # Check convergence here? (using diffRhoe)

        #
        Ehartree, Exc, Evtxc = update_from_rhoe!(Ham, Rhoe)

        descf = -sum( (Rhoe_in[:,1] .- Rhoe[:,1]).*(Vhartree + Vxc[:,1]) )*dVol

        Etot = Eband + deband + Ehartree + Exc + Ham.energies.NN + descf # entropy missed
    
        #
        diffEtot = abs(Etot - Etot_old)
        @printf("\nSCF: %5d %18.10f %10.5e %10.5e\n", iterSCF, Etot, diffEtot, diffRhoe)
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
    @printf("Eband    = %18.10f\n", Eband*2)
    @printf("deband   = %18.10f\n", deband*2)
    @printf("descf    = %18.10f\n", descf*2)
    @printf("-----------------------------\n")
    @printf("OneEle   = %18.10f\n", 2*(Eband + deband))
    @printf("Ehartree = %18.10f\n", 2*Ehartree)
    @printf("Exc      = %18.10f\n", 2*Exc)
    @printf("NN       = %18.10f\n", 2*Ham.energies.NN)
    @printf("-----------------------------\n")
    @printf("! Total  = %18.10f\n", 2*Etot)

    # TODO
    # Also print the Kohn-Sham orbital energies using similar format
    # as in pw.x


    if !is_converged
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    # Compare the energy using the usual formula (not using double-counting)
    energies = calc_energies(Ham, psiks)
    println("\nUsing original formula for total energy")
    println(energies)

    return
end


# Initialize Rhoe, potentials, reorthogonalize psiks with S
function _prepare_scf!(Ham, psiks)
    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham)
    # Update the potentials
    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, Rhoe, RhoeG )
    #
    # Reorthonormalize with S
    # FIXME: need to be included in rand_Blochwavefunc)
    # FIXME: We don't yet provide a function to do this
    #
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin - 1)*Nkpt
        ortho_sqrt_with_S!(Ham, psiks[ikspin])
    end
    return Ehartree, Exc, Evtxc
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

