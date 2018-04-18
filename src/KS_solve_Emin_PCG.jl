#
# Ham.potentials.V_Ps_loc should be initialized
# Ham.energies.NN should be calculated if needed
#
function KS_solve_Emin_PCG!( Ham::PWHamiltonian;
                             startingwfc=nothing, savewfc=true,
                             α_t=3e-5, NiterMax=200, verbose=true,
                             I_CG_BETA=2, ETOT_CONV_THR=1e-6, )

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    #
    # Initial wave function
    #
    if startingwfc == nothing
        srand(1234)
        psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
        psi = ortho_gram_schmidt(psi)
    else
        psi = startingwfc
    end

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    rhoe = calc_rhoe( pw, Focc, psi )
    update!(Ham, rhoe)

    #
    # Variables for PCG
    #
    d = zeros(Complex128, Ngwx, Nstates)
    g_old = zeros(Complex128, Ngwx, Nstates)
    d_old = zeros(Complex128, Ngwx, Nstates)
    Kg = zeros(Complex128, Ngwx, Nstates)
    Kg_old = zeros(Complex128, Ngwx, Nstates)
    β        = 0.0
    Etot_old = 0.0

    # Calculate energy at this psi
    energies = calc_energies(Ham, psi)
    Ham.energies = energies

    Etot     = energies.Total

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

        g = calc_grad( Ham, psi)
        Kg = Kprec(pw,g)

        if iter != 1
            if I_CG_BETA == 1
                β = real(sum(conj(g).*Kg))/real(sum(conj(g_old).*Kg_old))
            elseif I_CG_BETA == 2
                β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g_old).*Kg_old))
            elseif I_CG_BETA == 3
                β = real(sum(conj(g-g_old).*Kg))/real(sum(conj(g-g_old).*d))
            else
                β = real(sum(conj(g).*Kg))/real(sum((g-g_old).*conj(d_old)))
            end
        end
        if β < 0.0
            if verbose
                @printf("β is smaller than 0, setting it to zero\n")
            end
            β = 0.0
        end

        d = -Kprec(pw, g) + β * d_old

        psic = ortho_gram_schmidt(psi + α_t*d)
        rhoe = calc_rhoe( pw, Focc, psic )
        #
        update!(Ham, rhoe)

        gt = calc_grad(Ham, psic)

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        psi = psi[:,:] + α*d[:,:]

        # Update potentials
        psi = ortho_gram_schmidt(psi)
        rhoe = calc_rhoe( pw, Focc, psi )

        update!(Ham, rhoe)

        Ham.energies = calc_energies( Ham, psi )
        Etot = Ham.energies.Total
        diff = abs(Etot-Etot_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        end
        
        if diff < ETOT_CONV_THR
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
    psi = ortho_gram_schmidt(psi)
    Hr = psi' * op_H( Ham, psi )
    evals, evecs = eig(Hr)
    psi = psi*evecs

    if savewfc
        wfc_file = open("WFC.data","w")
        write(wfc_file,psi)
        close(wfc_file)
    end

    return

end
