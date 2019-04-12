function gen_Rhoe_aux( Ham:: Hamiltonian )
    return gen_Rhoe_aux( Ham.atoms, Ham.pw )
end

function gen_Rhoe_aux( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    # structure factor
    Sf = calc_strfact( atoms, pw )
    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    
    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    Rhoe_aux = zeros(Float64,Npoints)

    for ig = 1:Ng
        ip = idx_g2r[ig]
        for isp = 1:Nspecies
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] + Zvals[isp]*exp(-0.125*G2[ig]/eta^2)*Sf[ig,isp]/CellVolume
        end
    end
    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )
    return -Rhoe_aux*Npoints  # note the minus sign
end



"""
Solves Kohn-Sham problem using traditional self-consistent field (SCF)
iterations with density mixing.
"""
function KS_solve_SCF_03!( Ham::Hamiltonian ;
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
                           ETOT_CONV_THR=1e-6 )

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
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        Rhoe[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
    end

    if Nspin == 2
        @printf("\nInitial integ Rhoe up = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("\nInitial integ Rhoe dn = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("\nInitial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
    end

    Rhoe_aux = gen_Rhoe_aux( Ham )
    println("integ of Rhoe_aux = ", sum(Rhoe_aux)*dVol)

    V_aux_G = Poisson_solve( pw, Rhoe_aux )
    V_aux = real( G_to_R( pw, V_aux_G ) )
    println("integ of V_aux = ", sum(V_aux)*dVol)

    @assert Nspin == 1  # current limitation

    # update!(Ham, Rhoe) 
    # Update potentials here
    Rhoe_tot = Rhoe[:,1]
    Rhoe_neu = Rhoe_tot + Rhoe_aux
    println("integ Rhoe neu = ", sum(Rhoe_neu)*dVol)
    Ham.potentials.Hartree = real( G_to_R(pw, Poisson_solve(pw, Rhoe_neu)) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, Rhoe[:,1] )
    else  # VWN is the default
        Ham.potentials.XC[:,1] = calc_Vxc_VWN( Rhoe[:,1] )
    end
    # add potential due to Gaussian chgden
    #Ham.potentials.Hartree = Ham.potentials.Hartree #+ V_aux

    Ps_loc_orig = copy(Ham.potentials.Ps_loc)
    Ham.potentials.Ps_loc = Ham.potentials.Ps_loc - V_aux

    
    # Check integrated electrostatic potentials
    V_es_tot = Ham.potentials.Hartree + Ham.potentials.Ps_loc
    println("integ V es tot = ", sum(V_es_tot)*dVol)


    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ETHR_EVALS_LAST = 1e-6

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

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k


    for iter = 1:NiterMax

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

        if E_GAP_INFO && verbose && print_e_gap
            @printf("E gap = %18.10f\n", minimum(evals[idx_LUMO,:] - evals[idx_HOMO,:]))
        end

        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        Rhoe_new[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
        for ispin = 1:Nspin
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
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

        for rho in Rhoe
            if rho < eps()
                rho = 0.0
            end
        end

        # renormalize
        if check_rhoe
            integRhoe = sum(Rhoe)*dVol
            #@printf("After mixing: integRhoe = %18.10f\n", integRhoe)
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
            #@printf("After renormalize Rhoe: = %18.10f\n", integRhoe)
        end

        
        #update!( Ham, Rhoe )
        Ham.rhoe = Rhoe[:,:]
        Rhoe_tot = Rhoe[:,1]
        Rhoe_neu = Rhoe_tot + Rhoe_aux
        println("integ Rhoe neu = ", sum(Rhoe_neu)*dVol)
        Ham.potentials.Hartree = real( G_to_R(pw, Poisson_solve(pw, Rhoe_neu)) )
        if Ham.xcfunc == "PBE"
            Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, Rhoe[:,1] )
        else  # VWN is the default
            Ham.potentials.XC[:,1] = calc_Vxc_VWN( Rhoe[:,1] )
        end
        
        # add potential due to Gaussian chgden
        #Ham.potentials.Hartree = Ham.potentials.Hartree #+ V_aux


        #
        # Calculate energies
        #
        #Ham.energies = calc_energies( Ham, psiks )

        #
        # Kinetic energy
        #
        E_kin = 0.0
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            psi = psiks[ikspin]
            for ist = 1:Nstates
                psiKpsi = 0.0
                for igk = 1:Ngw[ik]
                    ig = idx_gw2g[ik][igk]
                    Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                    psiKpsi = psiKpsi + abs(psi[igk,ist])^2*Gw2
                end
                E_kin = E_kin + wk[ik]*Focc[ist,ikspin]*psiKpsi
            end
        end
        end
        E_kin = 0.5*E_kin

        #E_Hartree = 0.5*sum(Ham.potentials.Hartree.*Rhoe_neu)*dVol  # use G-space formula?
        E_Ps_loc  = sum(Ps_loc_orig.*Rhoe_tot)*dVol
        
        if Ham.xcfunc == "PBE"
            epsxc = calc_epsxc_PBE( Ham.pw, Rhoe )
        else
            epsxc = calc_epsxc_VWN( Rhoe )
        end
        E_xc = dot( epsxc, Rhoe_tot ) * dVol

        cRhoeG = conj(R_to_G(pw, Rhoe_neu)) #/sqrt(Npoints)
        V_HartreeG = R_to_G(pw, Ham.potentials.Hartree)
        #V_Ps_locG = R_to_G(pw, Ham.potentials.Ps_loc)

        #if Ham.xcfunc == "PBE"
        #    epsxc = calc_epsxc_PBE( Ham.pw, Ham.rhoe )
        #else
        #    epsxc = calc_epsxc_VWN( Ham.rhoe )
        #end
        #epsxcG = R_to_G(pw, epsxc)

        E_Hartree = 0.0
        #E_Ps_loc = 0.0
        for ig = 2:pw.gvec.Ng
            ip = pw.gvec.idx_g2r[ig]
            E_Hartree = E_Hartree + abs(cRhoeG[ip])^2/pw.gvec.G2[ig]
        #    E_Ps_loc = E_Ps_loc + real( V_Ps_locG[ip]*cRhoeG[ip] )
        end
        E_Hartree = 2*pi*CellVolume*E_Hartree/Npoints/Npoints/4.0
        #E_Hartree = 0.5*E_Hartree*dVol
        #E_Ps_loc = E_Ps_loc*dVol

        #E_xc = 0.0
        #for ig = 1:pw.gvec.Ng
        #    ip = pw.gvec.idx_g2r[ig]
        #    E_xc = E_xc + real( epsxcG[ip]*cRhoeG[ip] )
        #end
        #E_xc = E_xc*dVol


        if Ham.pspotNL.NbetaNL > 0
            E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
        else
            E_Ps_nloc = 0.0
        end

        Ham.energies.Kinetic = E_kin
        Ham.energies.Ps_loc  = E_Ps_loc
        Ham.energies.Ps_nloc = E_Ps_nloc
        Ham.energies.Hartree = E_Hartree
        Ham.energies.XC      = E_xc

        if use_smearing
            Ham.energies.mTS = Entropy
        end
        Etot = sum(Ham.energies)
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

        if diffE < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
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
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")
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
