function KS_solve_SCF!( Ham::HamiltonianGamma; kwargs... )
    KS_solve_SCF!( Ham, randn_BlochWavefunc(Ham); kwargs... )
    return
end

function KS_solve_SCF!(
    Ham::HamiltonianGamma, psis::BlochWavefuncGamma;
    NiterMax=150,
    betamix=0.2,
    startingrhoe=:random,
    verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    print_integ_rhoe=false,
    check_rhoe=false,
    use_smearing=false,
    mix_method="simple",
    mixdim=5,
    kT=1e-3,
    update_psi="LOBPCG",
    cheby_degree=8,
    etot_conv_thr=1e-6,
    ethr_evals_last=1e-5,
    starting_magn=nothing 
)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = 1
    Nkspin = Nspin*Nkpt
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    pspots = Ham.pspots
    electrons = Ham.electrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nelectrons = Ham.electrons.Nelectrons
    wk = [1.0]
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
        @printf("Density mixing with betamix = %10.5f\n", betamix)
        if use_smearing
            @printf("Smearing = %f\n", kT)
        end
        println("") # blank line before SCF iteration info
    end


    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_new = zeros(Float64,Npoints,Nspin)
    if startingrhoe == :gaussian
        error("starting rhoe gaussian is not yet implemented")
        if Nspin == 1
            Rhoe[:,1] = guess_rhoe( Ham )
        else
            Rhoe = guess_rhoe_atomic( Ham, starting_magn=starting_magn )
        end
    else
        calc_rhoe!( Ham, psis, Rhoe )
    end

    if Nspin == 2 && verbose
        @printf("Initial integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Initial integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Initial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
        println("")
    end

    if mix_method in ("anderson", "broyden")
        df = zeros(Float64,Npoints*Nspin, mixdim)
        dv = zeros(Float64,Npoints*Nspin, mixdim)
    
    elseif mix_method == "linear_adaptive"
        betav = betamix*ones(Float64, Npoints*Nspin)
        df = zeros(Float64, Npoints*Nspin)
    
    elseif mix_method in ("rpulay", "ppulay", "pulay")
        XX = zeros(Float64,Npoints*Nspin, mixdim)
        FF = zeros(Float64,Npoints*Nspin, mixdim)
        x_old = zeros(Float64,Npoints,Nspin)
        f_old = zeros(Float64,Npoints,Nspin)
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
        #if iterSCF == 1
        #    ethr = 0.1
        #elseif iterSCF == 2
        #    ethr = 0.01
        #else
        #    ethr = ethr/5.0
        #    ethr = max( ethr, ethr_evals_last )
        #end
        ethr = ethr_evals_last

        evals =
        diag_Emin_PCG!( Ham, psis, verbose=false, verbose_last=false, tol=ethr,
                        Nstates_conv=Nstates_occ )

        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        calc_rhoe!( Ham, psis, Rhoe_new )

        for ispin = 1:Nspin
            diffRhoe[ispin] = sum(abs.(Rhoe_new[:,ispin] - Rhoe[:,ispin]))/Npoints
        end

        # check norm of
        if check_rhoe
            integRhoe = sum(Rhoe_new)*dVol
            @printf("Before: integRhoe_new = %18.10f\n", integRhoe)
            Rhoe_new = Nelectrons/integRhoe * Rhoe_new
            integRhoe = sum(Rhoe_new)*dVol
            @printf("After renormalize Rhoe_new: = %18.10f\n", integRhoe)
        end

        if mix_method == "simple"

            Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        elseif mix_method == "linear_adaptive"

            mix_adaptive!( Rhoe, Rhoe_new, betamix, betav, df )

        elseif mix_method == "broyden"

            mix_broyden!( Rhoe, Rhoe_new, betamix, iterSCF, mixdim, df, dv )

        elseif mix_method == "pulay"
        
            mix_pulay!( Rhoe, Rhoe_new, betamix, XX, FF, iterSCF, mixdim, x_old, f_old )

        elseif mix_method == "rpulay"
            
            mix_rpulay!( Rhoe, Rhoe_new, betamix, XX, FF, iterSCF, mixdim, x_old, f_old )
            # result is in Rhoe

        elseif mix_method == "ppulay"
            
            #XXX We fix the period to be 3 here

            mix_ppulay!( Rhoe, Rhoe_new, betamix, XX, FF, iterSCF, mixdim, 3, x_old, f_old )

        
        elseif mix_method == "anderson"

            mix_anderson!( Rhoe, Rhoe_new, betamix, df, dv, iterSCF, mixdim )
        
        else
            error(@sprintf("Unknown mix_method = %s\n", mix_method))

        end

        for i in 1:length(Rhoe)
            if Rhoe[i] < eps()
                Rhoe[i] = eps()
            end
        end

        if Nspin == 2
            magn_den = Rhoe[:,1] - Rhoe[:,2]
        end

        # renormalize
        if check_rhoe
            integRhoe = sum(Rhoe)*dVol
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
        end

        update!( Ham, Rhoe )
        # Calculate energies
        Ham.energies = calc_energies( Ham, psis )
        if use_smearing
            Ham.energies.mTS = Entropy
        end
        Etot = sum(Ham.energies)
        diffE = abs( Etot - Etot_old )

        if verbose
            if Nspin == 1
                @printf("SCF: %5d %18.10f %12.5e %12.5e\n", iterSCF, Etot, diffE, diffRhoe[1] )
                if print_integ_rhoe
                    @printf("integ Rhoe = %18.10f\n", sum(Rhoe)*dVol)
                end
            else
                @printf("SCF: %5d %18.10f %12.5e [%12.5e,%12.5e]\n", iterSCF, Etot,
                    diffE, diffRhoe[1], diffRhoe[2] )
                if print_integ_rhoe
                    magn_den = Rhoe[:,1] - Rhoe[:,2]
                    @printf("integ Rhoe spin up = %18.10f\n", sum(Rhoe[:,1])*dVol) 
                    @printf("integ Rhoe spin dn = %18.10f\n", sum(Rhoe[:,2])*dVol) 
                    @printf("integ magn_den = %18.10f\n", sum(magn_den)*dVol) 
                end
            end     
        end

        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            if verbose
                @printf("\nSCF is converged in iter: %d\n", iterSCF)
            end
            break
        end
        #
        Etot_old = Etot

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

    return
end