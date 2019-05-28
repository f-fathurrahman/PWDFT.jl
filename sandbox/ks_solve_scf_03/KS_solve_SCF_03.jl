using SpecialFunctions: erf, besselj0
using QuadGK


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
                           etot_conv_thr=1e-6 )

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
        Rhoe[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin, renormalize=false )
    end

    if Nspin == 2
        @printf("\nInitial integ Rhoe up = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("\nInitial integ Rhoe dn = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("\nInitial integ magn_den = %18.10f\n", sum(Rhoe[:,1] - Rhoe[:,2])*dVol)
    end

    Rhoe_aux, E_self = gen_Rhoe_aux_G( Ham.atoms, Ham.pw, Ham.pspots )
    println("integ of Rhoe_aux = ", sum(Rhoe_aux)*dVol)
    println("E self = ", E_self)

    # Check integrated electrostatic potentials
    V_es_tot = Ham.potentials.Hartree + Ham.potentials.Ps_loc
    int_V_es = sum(V_es_tot)*dVol
    println("integ V es tot = ", sum(V_es_tot)*dVol)


    V_Ps_loc_short = gen_V_Ps_loc_short( Ham )

    V_aux_G = Poisson_solve( pw, Rhoe_aux )
    #V_aux_G[1] = 2*sqrt(1.0/pi)*Nelectrons*dVol
    ##V_aux_G[1] = Nelectrons
    V_aux = real( G_to_R( pw, V_aux_G ) )

    V_aux_v2 = gen_V_aux_G( Ham.atoms, Ham.pw, Ham.pspots )
    
    #V_aux_v2 = Ham.potentials.Ps_loc - V_Ps_loc_short

    println("Some V_aux: ")
    for i = 1:5
        dr = norm(Ham.pw.r[:,i])
        if dr <= eps()
            dr = eps()
        end
        V_erf = -Ham.atoms.Zvals[1]/dr*erf(dr)
        @printf("%5d %18.10f %18.10f %18.10f\n", i, V_aux[i], V_aux_v2[i], abs(V_erf-V_aux[i]) )
    end
    println("minimum V_aux = ", minimum(V_aux))
    
    Vzero = -2*sqrt(1.0/pi)*Ham.atoms.Zvals[1]
    println("-Zvals*erf(r)/r at zero = ", Vzero)  # for the first species

    println("diff Vzero = ", (Vzero - minimum(V_aux))*dVol)


    println("integ of V_Ps_loc_short = ", sum(V_Ps_loc_short)*dVol)
    println("integ of V_aux = ", sum(V_aux)*dVol)
    println("integ of V_Ps_loc = ", sum(Ham.potentials.Ps_loc)*dVol)

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

    Ps_loc_orig = copy(Ham.potentials.Ps_loc)
    #Ham.potentials.Ps_loc = Ham.potentials.Ps_loc - V_aux
    
    #Ham.potentials.Ps_loc = gen_V_Ps_loc_screened( Ham.atoms, Ham.pw, Ham.pspots )
    
    Ham.potentials.Ps_loc = V_Ps_loc_short
    #V_aux = Ps_loc_orig - V_Ps_loc_short

    #println("integ of V_Ps_loc = ", sum(Ham.potentials.Ps_loc)*dVol)

    #exit()

    
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

        Rhoe_new[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin, renormalize=false )
        for ispin = 1:Nspin
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
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
        
        @printf("integ Rhoe neu = %18.10f\n", sum(Rhoe_neu)*dVol)
        
        V_HartreeG = Poisson_solve(pw, Rhoe_neu)
        Ham.potentials.Hartree = real( G_to_R(pw, V_HartreeG) )

        # Check integrated electrostatic potentials
        V_es_tot = Ham.potentials.Hartree + Ham.potentials.Ps_loc
        println("integ V es tot = ", sum(V_es_tot)*dVol)

        if Ham.xcfunc == "PBE"
            Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, Rhoe[:,1] )
        else  # VWN is the default
            Ham.potentials.XC[:,1] = calc_Vxc_VWN( Rhoe[:,1] )
        end


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
        #E_Ps_loc  = sum(Ps_loc_orig[2:end].*Rhoe_tot[2:end])*dVol
        

        cRhoeG = conj(R_to_G(pw, Rhoe_neu))/Npoints
        cRhoeG_tot = conj(R_to_G(pw, Rhoe_tot))/Npoints
        V_Ps_locG = R_to_G(pw, Ps_loc_orig)

        if Ham.xcfunc == "PBE"
            epsxc = calc_epsxc_PBE( Ham.pw, Ham.rhoe )
        else
            epsxc = calc_epsxc_VWN( Ham.rhoe )
        end
        E_xc = dot( epsxc, Rhoe_tot ) * dVol

        E_Hartree = 0.0
        E_Ps_loc = 0.0
        for ig = 2:pw.gvec.Ng
            ip = pw.gvec.idx_g2r[ig]
            #E_Hartree = E_Hartree + 4*pi*abs(cRhoeG[ip])^2/pw.gvec.G2[ig]
            #E_Hartree = E_Hartree + real( (V_HartreeG[ip] - V_aux_G[ip]) * cRhoeG_tot[ip] )
            #E_Hartree = E_Hartree + 4*pi*abs(cRhoeG_tot[ip])^2 / pw.gvec.G2[ig]
            E_Ps_loc = E_Ps_loc + real( V_Ps_locG[ip] * cRhoeG_tot[ip] )
        end
        #E_Hartree = 0.5*( E_Hartree - V_aux_G[1] )*dVol #/CellVolume #- E_self
        #E_Hartree = 0.5*E_Hartree #*dVol
        E_Ps_loc = E_Ps_loc*dVol

        #E_Hartree = 0.5*sum( (Ham.potentials.Hartree - V_aux).*(Rhoe_neu - Rhoe_aux) )*dVol #- E_self
        #E_Hartree = 0.5*sum( (Ham.potentials.Hartree - V_aux).*Rhoe_tot )*dVol #- E_self
        
        E_Hartree = 0.5*sum( Ham.potentials.Hartree.*Rhoe_neu )*dVol #- E_self
        
        #E_Hartree = 0.5*sum( (Ham.potentials.Hartree - V_aux).*(Rhoe_neu - Rhoe_aux) )*dVol

        println("E_Hartree_neu = ", 0.5*sum( Ham.potentials.Hartree.*Rhoe_neu )*dVol)
        println("E_self        = ", E_self)
        println("E_Hartree     = ", E_Hartree)
        println("E_Hartree_old = ", 0.5*sum( (Ham.potentials.Hartree - V_aux).*(Rhoe_neu - Rhoe_aux) )*dVol)
        println("E_self num    = ", 0.5*sum( (V_aux).*(Rhoe_aux) )*dVol)

        println("E neu - aux   = ", 0.5*sum( Ham.potentials.Hartree.*Rhoe_neu )*dVol - sum( (V_aux).*(Rhoe) )*dVol)
        println("E neu-aux-num = ", -0.5*sum( (V_aux).*(Rhoe_aux) )*dVol + 0.5*sum( Ham.potentials.Hartree.*Rhoe_neu )*dVol - sum( (V_aux).*(Rhoe) )*dVol)

        println("E H Rhoe neu  = ", 0.5*sum( (Ham.potentials.Hartree).*(Rhoe_aux) )*dVol)
        
        E_alphat = calc_E_alphat(Ham.atoms, Ham.pw, Ham.pspots )
        
        println("E_alphat = ", E_alphat)

        E_Hartree = 0.5*sum( (Ham.potentials.Hartree - V_aux).*(Rhoe_neu - Rhoe_aux) )*dVol - E_alphat # override
        E_Ps_loc = sum( Ham.potentials.Ps_loc.*Rhoe )*dVol + sum( (V_aux).*(Rhoe) )*dVol
        #exit()

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

        if diffE < etot_conv_thr
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
