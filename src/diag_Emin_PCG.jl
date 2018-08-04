"""
Find eigenvalues and eigenvectors of `Ham::Hamiltonian` using `X0`
as initial guess for eigenvectors using preconditioned CG method.
    
**IMPORTANT** `X0` must be orthonormalized before.
"""
function diag_Emin_PCG( Ham::Hamiltonian, X0::Array{ComplexF64,2};
                        tol=1e-5, NiterMax=100, verbose=false,
                        verbose_last=false, Nstates_conv=0,
                        tol_ebands=1e-4,
                        α_t=3e-5, I_CG_BETA = 2 )

    ik = Ham.ik
    pw = Ham.pw
    Focc = Ham.electrons.Focc

    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(X0)[2]
    Ngw_ik = pw.gvecw.Ngw[ik]

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    psi = copy(X0)

    #
    # Variabls for PCG
    #
    d = zeros(ComplexF64, Ngw_ik, Nstates)
    g_old = zeros(ComplexF64, Ngw_ik, Nstates)
    d_old = zeros(ComplexF64, Ngw_ik, Nstates)
    Kg = zeros(ComplexF64, Ngw_ik, Nstates)
    Kg_old = zeros(ComplexF64, Ngw_ik, Nstates)
    psic = zeros(ComplexF64, Ngw_ik, Nstates)
    gt = zeros(ComplexF64, Ngw_ik, Nstates)
    
    β = 0.0
    
    Ebands_old = 0.0

    Hr = Hermitian( psi' * op_H( Ham, psi ) )
    
    evals = eigvals(Hr)
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false

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

        Hr = Hermitian( psi' * op_H( Ham, psi ) )
        
        evals = eigvals(Hr)
        Ebands = sum(evals)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        nconv = length( findall( devals .< tol ) )

        diff = abs(Ebands-Ebands_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diff)
            for ist = 1:Nstates
                @printf("evals[%3d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
            end
            @printf("iter %d nconv = %d\n", iter, nconv)
        end
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Convergence is achieved based on nconv\n")
            end
            break
        end
        if diff <= tol_ebands*Nstates
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Convergence is achieved based on tol_ebands*Nstates\n")
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands
    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_Emin_PCG is not converged after %d iterations\n", NiterMax)
    end

    psi = ortho_sqrt(psi)
    Hr = Hermitian( psi' * op_H(Ham, psi) )
    evals, evecs = eigen(Hr)
    # Sort
    idx_sorted = sortperm(evals)
    evals = evals[idx_sorted]
    psi = psi*evecs[:,idx_sorted]

    if verbose_last || verbose
        @printf("\nEigenvalues from diag_Emin_PCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end

    return evals, psi
end
