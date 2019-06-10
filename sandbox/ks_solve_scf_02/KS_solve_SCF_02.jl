using LinearAlgebra: norm

include("../guess_rhoe_atomic/guess_rhoe_atomic.jl")

function KS_solve_SCF_02!(
    Ham::Hamiltonian ;
    startingwfc=:random, savewfc=false,
    startingrhoe=:gaussian,
    betamix=0.2, NiterMax=100, verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    check_rhoe=false,
    use_smearing=false, kT=1e-3,
    update_psi="LOBPCG", cheby_degree=8,
    mix_method="simple", MIXDIM=5,
    print_e_gap=false,
    etot_conv_thr=1e-6
)

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    wk = Ham.pw.gvecw.kpoints.wk
    #
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    #
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints
    #
    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = electrons.Nspin
    #
    Nkspin = Nkpt*Nspin

    Nstates_occ = electrons.Nstates_occ

    #
    # Initial wave function
    #
    if startingwfc == :read
        psiks = read_psiks( Ham )
    else
        # generate random BlochWavefunc
        psiks = rand_BlochWavefunc( Ham )
    end

    E_GAP_INFO = false
    if Nstates_occ < Nstates
        E_GAP_INFO = true
        if Nspin == 2
            idx_HOMO = max(round(Int64,Nstates_occ/2),1)
            idx_LUMO = idx_HOMO + 1
        else
            idx_HOMO = Nstates_occ
            idx_LUMO = idx_HOMO + 1
        end
    end

    if Ham.sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)

    if Nspin == 2
        idx_HOMO = max(round(Int64,Nstates_occ/2),1)
        idx_LUMO = idx_HOMO + 1
        Focc[idx_HOMO,1:Nkpt] = Focc[idx_HOMO,1:Nkpt] .+ 0.5
        Focc[idx_HOMO,Nkpt+1:2*Nkpt] = Focc[idx_HOMO,Nkpt+1:2*Nkpt] .- 0.5
    end

    if startingrhoe == :gaussian && startingwfc == :random
        #@assert Nspin == 1
        #Rhoe[:,1] = guess_rhoe( Ham )
        Rhoe = guess_rhoe_atomic( Ham )
    else
        Rhoe[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
    end
    
    # Symmetrize Rhoe is needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
    end

    if Nspin == 2
        @printf("\nInitial integ Rhoe up = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("\nInitial integ Rhoe dn = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("\nInitial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
    end

    update!(Ham, Rhoe)

    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ETHR_EVALS_LAST = 1e-8

    ethr = 0.1
    
    if mix_method == "anderson"
        df = zeros(Float64,Npoints*Nspin,MIXDIM)
        dv = zeros(Float64,Npoints*Nspin,MIXDIM)
    
    elseif mix_method in ("rpulay", "rpulay_kerker", "ppulay", "pulay")
        XX = zeros(Float64,Npoints*Nspin,MIXDIM)
        FF = zeros(Float64,Npoints*Nspin,MIXDIM)
        x_old = zeros(Float64,Npoints*Nspin)
        f_old = zeros(Float64,Npoints*Nspin)
    end


    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("update_psi = %s\n", update_psi)
    @printf("\n")
    @printf("mix_method = %s\n", mix_method)
    if mix_method in ("rpulay", "rpulay_kerker", "anderson", "ppulay")
        @printf("MIXDIM = %d\n", MIXDIM)
    end
    @printf("Density mixing with betamix = %10.5f\n", betamix)
    if use_smearing
        @printf("Smearing = %f\n", kT)
    end
    println("") # blank line before SCF iteration info

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    CONVERGED = 0

    E_fermiSpin = zeros(Nspin)

    Etot = 0.0
    Etot_old = 999.0


    @assert Nspin == 1  # Current limitation

    Vin = zeros(Npoints)

    for iter = 1:NiterMax

        @printf("ethr = %18.10e\n", ethr)
        evals, ethr = do_diag!(Ham, psiks, update_psi, iter, Nstates_occ, ethr, ETHR_EVALS_LAST, cheby_degree)

        if E_GAP_INFO && verbose && print_e_gap
            @printf("E gap = %18.10f\n", minimum(evals[idx_LUMO,:] - evals[idx_HOMO,:]))
        end

        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        Rhoe_new[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
        # Symmetrize Rhoe is needed
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe_new )
        end


        for ispin = 1:Nspin
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
        end


        Eband = 0.0
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (Nspin - 1)*Nkpt
            for ist = 1:Nstates
                Eband = Eband + Focc[ist,ikspin]*wk[ik]*evals[ist,ikspin]
            end
        end
        end

        # 
        Vin[2:end] = Ham.potentials.Hartree[2:end] + Ham.potentials.XC[2:end,1]
        Vin[1] = Ham.potentials.XC[1,1]
        deband = -sum(Vin.*Rhoe_new[:,1])*dVol
        #deband = -sum(Vin.*Rhoe[:,1])*dVol

        if mix_method == "simple"
            for ispin = 1:Nspin
                Rhoe[:,ispin] = betamix*Rhoe_new[:,ispin] + (1-betamix)*Rhoe[:,ispin]
            end

        elseif mix_method == "simple_kerker"
            for ispin = 1:Nspin
                Rhoe[:,ispin] = Rhoe[:,ispin] + betamix*precKerker(pw, Rhoe_new[:,ispin] - Rhoe[:,ispin])
            end

        elseif mix_method == "pulay"
        
            Rhoe = reshape( mix_pulay!(
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )
            
            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end

        elseif mix_method == "rpulay"
        
            Rhoe = reshape( mix_rpulay!(
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )
            
            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end

        elseif mix_method == "ppulay"
        
            Rhoe = reshape( mix_ppulay!(
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, 3, x_old, f_old
                ), (Npoints,Nspin) )
            
            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end

        elseif mix_method == "rpulay_kerker"
        
            Rhoe = reshape( mix_rpulay_kerker!( pw,
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )
            
            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end
        
        elseif mix_method == "anderson"
            Rhoe[:,:] = mix_anderson!( Nspin, Rhoe, Rhoe_new, betamix, df, dv, iter, MIXDIM )
        
        else
            error(@sprintf("Unknown mix_method = %s\n", mix_method))

        end


        update!( Ham, Rhoe )

        descf = sum( (Rhoe[:,1] .- Rhoe_new[:,1]).*(Ham.potentials.Hartree .+ Ham.potentials.XC[:,1]) )*dVol

        EHartree = 0.5*sum(Ham.potentials.Hartree.*Rhoe[:,1])*dVol

        if Ham.xcfunc == "PBE"
            epsxc = calc_epsxc_PBE( Ham.pw, Ham.rhoe )
        else
            epsxc = calc_epsxc_VWN( Ham.rhoe )
        end
        Exc = sum(epsxc[:,1].*Rhoe[:,1])*dVol

        if use_smearing
            Ham.energies.mTS = Entropy
        end
        
        #Etot = sum(Ham.energies)
        Etot = Eband + deband + EHartree + Exc + Ham.energies.NN + descf + Ham.energies.PspCore # entropy missed

        @printf("Eband    = %18.10f\n", Eband)
        @printf("deband   = %18.10f\n", deband)
        @printf("OneEle   = %18.10f\n", Eband + deband)
        @printf("descf    = %18.10f\n", descf)
        @printf("EHartree = %18.10f\n", EHartree)
        @printf("Exc      = %18.10f\n", Exc)
        @printf("NN       = %18.10f\n", Ham.energies.NN)
        @printf("PspCore  = %18.10f\n", Ham.energies.PspCore)

        diffE = abs( Etot - Etot_old )

        if verbose
            if Nspin == 1
                @printf("\nSCF: %8d %18.10f %18.10e %18.10e\n",
                        iter, Etot, diffE, diffRhoe[1] )
                @printf("integ Rhoe = %18.10f\n", sum(Rhoe)*dVol)
            else
                @printf("SCF: %8d %18.10f %18.10e %18.10e %18.10e\n",
                        iter, Etot, diffE, diffRhoe[1], diffRhoe[2] )
                magn_den = Rhoe[:,1] - Rhoe[:,2]
                @printf("integ Rhoe spin up = %18.10f\n", sum(Rhoe[:,1])*dVol) 
                @printf("integ Rhoe spin dn = %18.10f\n", sum(Rhoe[:,2])*dVol) 
                @printf("integ magn_den = %18.10f\n", sum(magn_den)*dVol) 
            end
        
        end

        if diffE < etot_conv_thr
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
                println("")
                @printf("Eband    = %18.10f\n", Eband)
                @printf("deband   = %18.10f\n", deband)
                @printf("OneEle   = %18.10f\n", Eband + deband)
                @printf("descf    = %18.10f\n", descf)
                @printf("EHartree = %18.10f\n", EHartree)
                @printf("Exc      = %18.10f\n", Exc)
                @printf("NN       = %18.10f\n", Ham.energies.NN)
            end
            break
        end
        #
        Etot_old = Etot

        flush(stdout)
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




function do_diag!(
    Ham,
    psiks,
    update_psi,
    iter,
    Nstates_occ,
    ethr, ETHR_EVALS_LAST,
    cheby_degree
)

    # determine convergence criteria for diagonalization
    if iter == 1
        ethr = 0.1
    elseif iter == 2
        ethr = 0.01
    else
        ethr = ethr/5.0
        ethr = max( ethr, ETHR_EVALS_LAST )
    end

    if update_psi == "LOBPCG"
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false,
                      Nstates_conv=Nstates_occ )

    elseif update_psi == "davidson"
        evals =
        diag_davidson!( Ham, psiks, verbose=false, verbose_last=false,
                        Nstates_conv=Nstates_occ )                

    elseif update_psi == "PCG"
        evals =
        diag_Emin_PCG!( Ham, psiks, verbose=false, verbose_last=false,
                        Nstates_conv=Nstates_occ )

    elseif update_psi == "CheFSI"
        evals =
        diag_CheFSI!( Ham, psiks, cheby_degree )

    else
        error( @sprintf("Unknown method for update_psi = %s\n", update_psi) )
    end

    return evals, ethr

end

