function alt1_KS_solve_SCF_spinpol!( Ham::Hamiltonian ;
                             startingwfc=nothing,
                             betamix = 0.7, NiterMax=100, verbose=false,
                             check_rhoe_after_mix=true,
                             update_psi="LOBPCG", cheby_degree=8,
                             ETOT_CONV_THR=1e-6 )
    
    @assert( Ham.pw.gvecw.kpoints.Nkpt == 1 )
    @assert( Ham.electrons.Nspin == 2 )

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Ns = pw.Ns
    Npoints = prod(Ns)
    dVol = pw.CellVolume/Npoints
    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ

    Nspin = 2  # HARDCODED !!
    Nkpt = 1 # HARDCODED

    # Random guess of wave function
    if startingwfc==nothing
        psiks = gen_rand_wavefunc(pw, electrons)
    else
        psiks = startingwfc
    end

    # Calculate electron density from this wave function and update Hamiltonian
    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
    end
    update!(Ham, Rhoe)

    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    magn_den = Rhoe[:,1] - Rhoe[:,2]

    println("integ total Rhoe = ", sum(Rhoe)*dVol)
    println("integ magn_den   = ", sum(magn_den)*dVol)

    evals = zeros(Float64,Nstates,Nspin)

    ETHR_EVALS_LAST = 1e-6

    ethr = 0.1

    MIXDIM = 5
    XX = zeros(Float64,Npoints*Nspin,MIXDIM)
    FF = zeros(Float64,Npoints*Nspin,MIXDIM)

    x_old = zeros(Float64,Npoints*Nspin)
    f_old = zeros(Float64,Npoints*Nspin)

    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("\n")
    @printf("Density mixing with beta = %10.5f\n", betamix)
    @printf("\n")

    CONVERGED = 0

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
            for ispin = 1:Nspin
                evals[:,ispin], psiks[ispin] =
                diag_LOBPCG( Ham, psiks[ispin], tol=ethr, verbose=false, verbose_last=false,
                             Nstates_conv = Nstates_occ )
            end
        #elseif update_psi == "davidson"
        #    evals, psi =
        #    diag_davidson( Ham, psi, tol=ethr, verbose=false, verbose_last=false,
        #                   Nstates_conv = Nstates_occ )
#
        #elseif update_psi == "PCG"
        #    evals, psi =
        #    diag_Emin_PCG( Ham, psi, tol=ethr, verbose=false, verbose_last=false,
        #                   Nstates_conv = Nstates_occ )
#
        #elseif update_psi == "CheFSI"
        #    
        #    ub, lb = get_ub_lb_lanczos( Ham, Nstates*2 )
        #    psi = chebyfilt( Ham, psi, cheby_degree, lb, ub)
        #    psi = ortho_sqrt( psi )

        else

            @printf("ERROR: Unknown method for update_psi = %s\n", update_psi)
            exit()
        end

        for ispin = 1:Nspin
            idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
            Rhoe_new[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
        end

        Rhoe = reshape( mix_rpulay!(
            reshape(Rhoe,(Npoints*Nspin)),
            reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
            ), (Npoints,Nspin) )
        
        magn_den = Rhoe[:,1] - Rhoe[:,2]
        println("integ total Rhoe = ", sum(Rhoe)*dVol)
        println("integ magn_den   = ", sum(magn_den)*dVol)

        for rho in Rhoe
            if rho < eps()
                rho = 0.0
            end
        end

        if check_rhoe_after_mix
            integRhoe = sum(Rhoe)*dVol
            #@printf("After mixing: integRhoe = %18.10f\n", integRhoe)
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
            #@printf("After renormalize Rhoe: = %18.10f\n", integRhoe)
        end

        update!( Ham, Rhoe )

        # Calculate energies
        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = abs( Etot - Etot_old )

        @printf("SCF: %8d %18.10f %18.10e\n", iter, Etot, diffE )

        if diffE < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        #
        Etot_old = Etot
    end

    # Eigenvalues are not calculated if using CheFSI.
    # We calculate them here.
    #if update_psi == "CheFSI"
    #    Hr = psi' * op_H( Ham, psi )
    #    evals = real(eigvals(Hr))
    #end

    Ham.electrons.ebands[:,:] = evals

    return

end
