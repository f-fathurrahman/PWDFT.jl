#
# Ham.potentials.V_Ps_loc should be initialized
# Ham.energies.NN should be calculated if needed
#
function KS_solve_Emin_PCG!( Ham::PWHamiltonian;
                             startingwfc=nothing, savewfc=true,
                             α_t=3e-5, NiterMax=200, verbose=true,
                             I_CG_BETA=2, ETOT_CONV_THR=1e-6 )

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt

    psik = Array{Array{Complex128,2},1}(Nkpt)

    #
    # Initial wave function
    #
    if startingwfc == nothing
        srand(1234)
        for ik = 1:Nkpt
            psi = rand(Complex128,Ngw[ik],Nstates)
            psik[ik] = ortho_gram_schmidt(psi)
        end
    else
        psi = startingwfc
    end

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    update!(Ham, rhoe)

    #
    # Variables for PCG
    #
    g = Array{Array{Complex128,2},1}(Nkpt)
    d = Array{Array{Complex128,2},1}(Nkpt)
    g_old = Array{Array{Complex128,2},1}(Nkpt)
    d_old = Array{Array{Complex128,2},1}(Nkpt)
    Kg = Array{Array{Complex128,2},1}(Nkpt)
    Kg_old = Array{Array{Complex128,2},1}(Nkpt)
    psic = Array{Array{Complex128,2},1}(Nkpt)
    gt = Array{Array{Complex128,2},1}(Nkpt)
    #
    for ik = 1:Nkpt
        g[ik] = zeros(Complex128, Ngw[ik], Nstates)
        d[ik] = zeros(Complex128, Ngw[ik], Nstates)
        g_old[ik] = zeros(Complex128, Ngw[ik], Nstates)
        d_old[ik] = zeros(Complex128, Ngw[ik], Nstates)
        Kg[ik] = zeros(Complex128, Ngw[ik], Nstates)
        Kg_old[ik] = zeros(Complex128, Ngw[ik], Nstates)
        psic[ik] = zeros(Complex128, Ngw[ik], Nstates)
        gt[ik] = zeros(Complex128, Ngw[ik], Nstates)
    end
    
    β = zeros(Nkpt)
    α = zeros(Nkpt)

    Etot_old = 0.0

    # Calculate energy at this psi
    energies = calc_energies(Ham, psik)
    Ham.energies = energies

    Etot     = energies.Total

    CONVERGED = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG\n")
        @printf("-------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_trh  = %e\n", ETOT_CONV_THR)
        if I_CG_BETA == 1
            @printf("Using Fletcher-Reeves formula for CG_BETA\n")
        elseif I_CG_BETA == 2
            @printf("Using Polak-Ribiere formula for CG_BETA\n")
        elseif I_CG_BETA == 3
            @printf("Using Hestenes-Stiefeld formula for CG_BETA\n")
        else
            @printf("Using Dai-Yuan formula for CG_BETA\n")
        end
        @printf("\n")
    end


    for iter = 1:NiterMax

        for ik = 1:Nkpt

            Ham.ik = ik
            g[ik] = calc_grad( Ham, psik[ik] )
            Kg[ik] = Kprec( Ham.ik, pw, g[ik] )

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                if I_CG_BETA == 1
                    β[ik] =
                    real(sum(conj(g[ik]).*Kg[ik]))/real(sum(conj(g_old[ik]).*Kg_old[ik]))
                elseif I_CG_BETA == 2
                    β[ik] =
                    real(sum(conj(g[ik]-g_old[ik]).*Kg[ik]))/real(sum(conj(g_old[ik]).*Kg_old[ik]))
                elseif I_CG_BETA == 3
                    β[ik] =
                    real(sum(conj(g[ik]-g_old[ik]).*Kg[ik]))/real(sum(conj(g[ik]-g_old[ik]).*d[ik]))
                else
                    β[ik] =
                    real(sum(conj(g[ik]).*Kg[ik]))/real(sum((g[ik]-g_old[ik]).*conj(d_old[ik])))
                end
            end
            if β[ik] < 0.0
                if verbose
                    @printf("β is smaller than 0, setting it to zero\n")
                end
                β[ik] = 0.0
            end

            d[ik] = -Kg[ik] + β[ik] * d_old[ik]

            psic[ik] = ortho_gram_schmidt(psik[ik] + α_t*d[ik])
        end

        rhoe = calc_rhoe( pw, Focc, psic )
        #
        update!(Ham, rhoe)

        for ik = 1:Nkpt
            Ham.ik = ik
            gt[ik] = calc_grad(Ham, psic[ik])

            denum = real(sum(conj(g[ik]-gt[ik]).*d[ik]))
            if denum != 0.0
                α[ik] = abs( α_t*real(sum(conj(g[ik]).*d[ik]))/denum )
            else
                α[ik] = 0.0
            end

            # Update wavefunction
            psik[ik] = psik[ik] + α[ik]*d[ik]

            # Update potentials
            psik[ik] = ortho_gram_schmidt(psik[ik])
        end

        rhoe = calc_rhoe( pw, Focc, psik )

        update!(Ham, rhoe)

        Ham.energies = calc_energies( Ham, psik )
        Etot = Ham.energies.Total
        diff = abs(Etot-Etot_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        end
        
        if diff < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("CONVERGENCE ACHIEVED\n")
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot
    end

    # Calculate eigenvalues
    for ik = 1:Nkpt
        Ham.ik = ik
        psik[ik] = ortho_gram_schmidt(psik[ik])
        Hr = psik[ik]' * op_H( Ham, psik[ik] )
        evals, evecs = eig(Hr)
        Ham.electrons.ebands[:,ik] = real(evals[:])
        psik[ik] = psik[ik]*evecs
    end

    if savewfc
        for ik = 1:Nkpt
            wfc_file = open("WFC_k_"*string(ik)*".data","w")
            write( wfc_file, psik[ik] )
            close( wfc_file )
        end
    end

    return

end
