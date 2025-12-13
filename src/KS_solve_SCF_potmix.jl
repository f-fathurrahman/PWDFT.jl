function KS_solve_SCF_potmix!( Ham::Hamiltonian; kwargs... )
    KS_solve_SCF_potmix!( Ham, rand_BlochWavefunc(Ham); kwargs... )
    return
end

"""
    KS_solve_SCF_potmix!( Ham, kwargs... )

Solve Kohn-Sham problem using traditional self-consistent field (SCF)
iterations with potential mixing. Hartree and XC potentials are
mixed separately.

Most arguments in `KS_solve_SCF!` are supported however not all mixing
methods are implemented.
"""
function KS_solve_SCF_potmix!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    NiterMax::Int64=150,
    betamix::Float64=0.2,
    startingrhoe::Symbol=:gaussian,
    verbose::Bool=true,
    print_final_ebands::Bool=false,
    print_final_energies::Bool=true,
    savewfc::Bool=false,
    use_smearing::Bool=false,
    mix_method::String="simple",
    mixdim::Int64=5,
    kT::Float64=1e-3,
    update_psi::String="LOBPCG",
    cheby_degree::Int64=8,
    etot_conv_thr::Float64=1e-6,
    ethr_evals_last::Float64=1e-5,
    starting_magnetization::Union{Vector{Float64},Nothing}=nothing 
)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nspin*Nkpt
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    pspots = Ham.pspots
    electrons = Ham.electrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nelectrons = Ham.electrons.Nelectrons
    wk = Ham.pw.gvecw.kpoints.wk
    Nstates_occ = electrons.Nstates_occ
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    if verbose
        @printf("\n")
        @printf("Self-consistent iteration begins ...\n")
        @printf("update_psi = %s\n", update_psi)
        @printf("mix_method = %s\n", mix_method)
        if mix_method in ("rpulay", "anderson", "ppulay", "broyden")
            @printf("mixdim = %d\n", mixdim)
        end
        @printf("Potential mixing with betamix = %10.5f\n", betamix)
        if use_smearing
            @printf("Smearing = %f\n", kT)
        end
        println("") # blank line before SCF iteration info
    end


    Rhoe = zeros(Float64,Npoints,Nspin)
    if startingrhoe == :gaussian
        if Nspin == 1
            Rhoe[:,1] = guess_rhoe( Ham )
        else
            Rhoe = guess_rhoe_atomic( Ham, starting_magnetization=starting_magnetization )
        end
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end

    if Nspin == 2 && verbose
        @printf("Initial integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Initial integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Initial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
        println("")
    end

    Vxc_inp = zeros(Float64, Npoints, Nspin)
    VHa_inp = zeros(Float64, Npoints)

    if mix_method == "broyden"
        df_VHa = zeros(Float64, Npoints, mixdim)
        df_Vxc = zeros(Float64, Npoints*Nspin, mixdim)
        #
        dv_VHa = zeros(Float64, Npoints, mixdim)
        dv_Vxc = zeros(Float64, Npoints*Nspin, mixdim)

    elseif mix_method == "linear_adaptive"
        betav_Vxc = betamix*ones(Float64, Npoints*Nspin)
        df_Vxc = zeros(Float64, Npoints*Nspin)
        #
        betav_VHa = betamix*ones(Float64, Npoints)
        df_VHa = zeros(Float64, Npoints)
    end

    update!(Ham, Rhoe)

    Ham.energies.NN = calc_E_NN(atoms)

    evals = zeros(Nstates,Nkspin)

    Etot_old = 0.0

    Nconverges = 0

    ethr = 0.1
    diffRhoe = ones(Nspin)
    diffPot = ones(Nspin)
    Rhoe_old = zeros(Float64,Npoints,Nspin)
    E_fermi = 0.0

    if verbose
        if Nspin == 1
            @printf("--------------------------------------------------------------\n")
            @printf("              iter            E            ΔE           Δρ\n")
            @printf("--------------------------------------------------------------\n")
        else
            @printf("----------------------------------------------------------------------------\n")
            @printf("              iter            E            ΔE                  Δρ\n")
            @printf("----------------------------------------------------------------------------\n")
        end
    end

    for iterSCF = 1:NiterMax

        # determine convergence criteria for diagonalization
        if iterSCF == 1
            ethr = 0.1
        elseif iterSCF == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end


        if update_psi == "LOBPCG"

            evals =
            diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                          Nstates_conv=Nstates_occ )

        elseif update_psi == "davidson"

            evals =
            diag_davidson!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                            Nstates_conv=Nstates_occ )                

        elseif update_psi == "PCG"

            evals =
            diag_Emin_PCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                            Nstates_conv=Nstates_occ )

        elseif update_psi == "CheFSI"

            evals =
            diag_CheFSI!( Ham, psiks, cheby_degree )

        else
            error( @sprintf("Unknown method for update_psi = %s\n", update_psi) )
        end


        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Ham.energies.mTS = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.E_fermi = E_fermi
            Ham.electrons.Focc = copy(Focc)
        end

        calc_rhoe!( Ham, psiks, Rhoe )
        # Ensure that there is no negative rhoe
        for i in 1:Npoints*Nspin
            if Rhoe[i] < eps()
                Rhoe[i] = eps()
            end
        end

        # Save the old (input) potential
        Vxc_inp[:,:] = Ham.potentials.XC
        VHa_inp[:] = Ham.potentials.Hartree

        # Update potentials
        update!(Ham, Rhoe)

        # Now Ham.potentials contains new (output) potential
        
        # Calculate energies
        Ham.energies = calc_energies(Ham, psiks)
        if use_smearing
            Ham.energies.mTS = Entropy
        end
        Etot = sum(Ham.energies)

        # Calculate MAE
        for ispin = 1:Nspin
            ss = 0.0
            for ip in 1:Npoints
                ss = ss + abs(Rhoe[ip,ispin] - Rhoe_old[ip,ispin])
            end
            diffRhoe[ispin] = ss/Npoints
        end

        diffEtot = abs(Etot - Etot_old)

        if verbose
            if Nspin == 1
                @printf("SCF_potmix: %5d %18.10f %12.5e %12.5e\n", iterSCF, Etot, diffEtot, diffRhoe[1])
            else
                @printf("SCF_potmix: %5d %18.10f %12.5e [%12.5e,%12.5e]\n",
                    iterSCF, Etot, diffEtot, diffRhoe[1], diffRhoe[2])
            end
        end

        if diffEtot < etot_conv_thr
            Nconverges = Nconverges + 1
        else  # reset Nconverges
            Nconverges = 0
        end

        if (Nconverges >= 2)
            if verbose
                @printf("\nSCF_potmix is converged in iter: %d\n", iterSCF)
            end
            break
        end
        
        Etot_old = Etot
        Rhoe_old = copy(Rhoe)

        # Mix potentials (Hartree and XC, separately)

        if mix_method == "broyden"
            mix_broyden!( Ham.potentials.Hartree, VHa_inp, betamix, iterSCF, mixdim, df_VHa, dv_VHa )
            mix_broyden!( Ham.potentials.XC, Vxc_inp, betamix, iterSCF, mixdim, df_Vxc, dv_Vxc )

        elseif mix_method == "linear_adaptive"
            mix_adaptive!( Ham.potentials.Hartree, VHa_inp, betamix, betav_VHa, df_VHa )
            mix_adaptive!( Ham.potentials.XC, Vxc_inp, betamix, betav_Vxc, df_Vxc )

        else
            # simple mixing
            Ham.potentials.Hartree = betamix*Ham.potentials.Hartree + (1-betamix)*VHa_inp
            Ham.potentials.XC = betamix*Ham.potentials.XC + (1-betamix)*Vxc_inp
        end


        # Don't forget to update the total local potential
        for ispin = 1:Nspin
            for ip = 1:Npoints
                Ham.potentials.Total[ip,ispin] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                                 Ham.potentials.XC[ip,ispin]
            end
        end

        flush(stdout)
    end

    if Nconverges < 2
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    Ham.electrons.ebands = evals

    if use_smearing && verbose
        @printf("\nFermi energy = %18.10f Ha = %18.10f eV\n", E_fermi, E_fermi*2*Ry2eV)
    end

    if Nspin == 2 && verbose
        @printf("\n")
        @printf("Final integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Final integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Final integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
    end

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