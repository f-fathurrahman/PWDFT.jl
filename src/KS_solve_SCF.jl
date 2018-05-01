function KS_solve_SCF!( Ham::PWHamiltonian ;
                       startingwfc=nothing, savewfc=true,
                       β = 0.5, NiterMax=100, verbose=false,
                       check_rhoe_after_mix=false,
                       update_psi="LOBPCG", cheby_degree=8,
                       ETOT_CONV_THR=1e-6 )

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates

    psik = Array{Array{Complex128,2},1}(Nkpt)

    #
    # Random guess of wave function
    #
    if startingwfc==nothing
        srand(1234)
        for ik = 1:Nkpt
            psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
            psik[ik] = ortho_gram_schmidt(psi)
        end
    else
        psik = startingwfc
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    update!(Ham, rhoe)

    Etot_old = 0.0
    rhoe_old = copy(rhoe)

    evals = zeros(Float64,Nstates,Nkpt)

    const ETHR_EVALS_LAST = 1e-6

    ethr = 0.1

    #
    # For mixing
    #
    MIXDIM = 4
    df = zeros(Float64,Npoints,MIXDIM)
    dv = zeros(Float64,Npoints,MIXDIM)

    for iter = 1:NiterMax

        if update_psi == "LOBPCG"
            for ik = 1:Nkpt
                Ham.ik = ik
                evals[:,ik], psik[ik] = diag_lobpcg( Ham, psik[ik], verbose_last=true )
            end

        elseif update_psi == "PCG"
            
            # determined convergence criteria for diagonalization
            if iter == 1
                ethr = 0.1
            elseif iter == 2
                ethr = 0.01
            else
                ethr = ethr/5.0
                ethr = max( ethr, ETHR_EVALS_LAST )
            end

            for ik = 1:Nkpt
                Ham.ik = ik
                evals[:,ik], psik[ik] = diag_Emin_PCG( Ham, psik[ik], TOL_EBANDS=ethr )
            end

        elseif update_psi == "CheFSI"
            
            for ik = 1:Nkpt
                Ham.ik = ik
                ub, lb = get_ub_lb_lanczos( Ham, Nstates*2 )
                psik[ik] = chebyfilt( Ham, psik[ik], cheby_degree, lb, ub)
                psik[ik] = ortho_gram_schmidt( psik[ik] )
            end

        end

        #
        rhoe_new = calc_rhoe( pw, Focc, psik )
        diffRho = norm(rhoe_new - rhoe)/Npoints

        #rhoe = β*rhoe_new[:] + (1-β)*rhoe[:]

        rhoe = andersonmix!( rhoe, rhoe_new, β, df, dv, iter, MIXDIM )

        for ip = 1:Npoints
            if rhoe[ip] < 1e-12
                rhoe[ip] = 1e-12
            end
        end

        if check_rhoe_after_mix
            integRhoe = sum(rhoe)*ΔV
            @printf("After mixing: integRho = %18.10f\n", integRhoe)
        end

        update!( Ham, rhoe )

        # Calculate energies
        Ham.energies = calc_energies( Ham, psik )
        Etot = Ham.energies.Total
        diffE = abs( Etot - Etot_old )

        #
        @printf("SCF: %8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )

        if diffE < ETOT_CONV_THR
            @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end
        #
        Etot_old = Etot
    end

    # Eigenvalues are not calculated if using CheFSI.
    # We calculate them here.
    if update_psi == "CheFSI"
        for ik = 1:Nkpt
            Ham.ik = ik
            Hr = psik[ik]' * op_H( Ham, psik[ik] )
            evals[:,ik] = real(eigvals(Hr))
        end
    end

    Ham.electrons.ebands = evals

    if savewfc
        for ik = 1:Nkpt
            wfc_file = open("WFC_k_"*string(ik)*".data","w")
            write( wfc_file, psik[ik] )
            close( wfc_file )
        end
    end

    return

end
