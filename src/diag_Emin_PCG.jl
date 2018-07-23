function diag_Emin_PCG( Ham::Hamiltonian, X0::Array{ComplexF64,2};
                        α_t = 3e-5, NiterMax=200, verbose=false,
                        verbose_last = false,
                        I_CG_BETA = 2, TOL_EBANDS=1e-6 )


    ik = Ham.ik
    pw = Ham.pw
    Focc = Ham.electrons.Focc

    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(X0)[2]
    Ngw_ik = pw.gvecw.Ngw[ik]

    psi = copy(X0)  # X0 is assumed to be already orthonormalized

    #
    # Variabls for PCG
    #
    d = zeros(ComplexF64, Ngw_ik, Nstates)
    g_old = zeros(ComplexF64, Ngw_ik, Nstates)
    d_old = zeros(ComplexF64, Ngw_ik, Nstates)
    Kg = zeros(ComplexF64, Ngw_ik, Nstates)
    Kg_old = zeros(ComplexF64, Ngw_ik, Nstates)
    β        = 0.0
    Ebands_old = 0.0

    Hr = psi' * op_H( Ham, psi )
    Ebands = sum( eigvals(Hermitian(Hr)) )

    for iter = 1:NiterMax

        g = calc_grad_evals( Ham, psi)
        Kg = Kprec( Ham.ik, pw, g )

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
            β = 0.0
        end

        d = -Kg + β * d_old

        psic = ortho_sqrt( psi + α_t*d )
        gt = calc_grad_evals( Ham, psic )

        denum = real(sum(conj(g-gt).*d))
        if denum != 0.0
            α = abs( α_t*real(sum(conj(g).*d))/denum )
        else
            α = 0.0
        end

        # Update wavefunction
        psi = psi[:,:] + α*d[:,:]

        # Update potentials
        psi = ortho_sqrt(psi)

        Hr = psi' * op_H( Ham, psi )
        Ebands = sum( eigvals(Hermitian(Hr)) )

        diff = abs(Ebands-Ebands_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diff)
        end
        if diff < TOL_EBANDS
            if verbose
                @printf("CONVERGENCE ACHIEVED\n")
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands
    end

    psi = ortho_sqrt(psi)
    Hr = Hermitian( psi' * op_H(Ham, psi) )
    evals, evecs = eigen(Hr)
    # Sort
    idx_sorted = sortperm(evals)
    evals = evals[idx_sorted]
    psi = psi*evecs[:,idx_sorted]

    if verbose_last
        for j = 1:Nstates
            @printf("eigval[%2d] = %18.10f\n", j, evals[j] )
        end
    end

    return evals, psi
end
