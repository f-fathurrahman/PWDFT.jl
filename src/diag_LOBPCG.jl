"""
Find eigenvalues and eigenvectors of `Ham::Hamiltonian` using `X0`
as initial guess for eigenvectors using
locally-optimal block preconditioned conjugate gradient method
of Knyazev.

**IMPORTANT** `X0` must be orthonormalized before.
"""
function diag_LOBPCG( Ham::Hamiltonian, X0::Array{ComplexF64,2};
                      tol=1e-5, NiterMax=100, verbose=false,
                      verbose_last=false, Nstates_conv=0 )

    pw = Ham.pw
    # get size info
    Nbasis = size(X0)[1]
    Nstates = size(X0)[2]

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    X = copy(X0)

    HX = op_H( Ham, X )

    nconv = 0
    iter = 1
    resnrm = ones(Nstates,1)

    P = zeros(ComplexF64,Nbasis,Nstates)
    HP = zeros(ComplexF64,Nbasis,Nstates)

    evals = zeros(Float64,Nstates)
    devals = ones(Float64,Nstates)
    evals_old = copy(evals)

    tfudge = 1e5

    IS_CONVERGED = false
    
    for iter = 1:NiterMax

        # Rayleigh quotient
        XHX = Hermitian(X'*HX)
        evals = eigvals(XHX)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        # Residuals
        R = HX - X*XHX
        for j = 1:Nstates
            resnrm[j] = norm( R[:,j] )
        end

        nconv = length( findall( devals .< tol ) )
        
        ilock = findall( devals .<= tol/tfudge )
        iact  = findall( devals .> tol/tfudge )

        Nlock = length(ilock)
        Nact = length(iact)

        if verbose
            @printf("LOBPCG iter = %8d, nconv=%d\n", iter, nconv)
            for ist = 1:Nstates
                @printf("evals[%d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
            end
            println("Nlock = ", Nlock, " Nact = ", Nact)
        end

        # Check for convergence
        if nconv >= Nstates_conv
            IS_CONVERGED = true
            if verbose || verbose_last
                println("LOBPCG convergence: nconv in iter ", iter)
            end
            break
        end

        # apply preconditioner
        W = Kprec( Ham.ik, pw, R )
        
        #
        # nlock == 0
        #
        HW = op_H( Ham, W )
        #
        C  = W'*W
        #C = 0.5*( C + C' )
        #
        R  = (cholesky(Hermitian(C))).U
        W  = W/R
        HW = HW/R
        #
        Q  = [X W]
        HQ = [HX HW]
        if iter > 1
            Q  = [Q P]
            HQ = [HQ HP]
        end

        T = Q'*(HQ)
        #T = 0.5*(T+T')
        T = Hermitian(T)
        
        G = Q'*Q
        #G = 0.5*(G+G')
        G = Hermitian(G)

        sd, S = eigen( T, G ) # evals, evecs
        U = S[:,1:Nstates]
        X = Q*U
        HX = HQ*U
        if iter > 1
            set2 = Nstates+1:2*Nstates
            set3 = 2*Nstates+1:3*Nstates
            P  = W*U[set2,:] + P*U[set3,:]
            HP = HW*U[set2,:] + HP*U[set3,:]
            C = P'*P
            #C = 0.5*(C + C')
            R = (cholesky(Hermitian(C))).U
            P = P/R
            HP = HP/R
        else
            P  = copy(W)
            HP = copy(HW)
        end

        iter = iter + 1
    end

    if !IS_CONVERGED
        @printf("\nWARNING: LOBPCG is not converged after %d iterations\n", NiterMax)
    end

    S = X'*HX
    S = Hermitian(0.5*(S+S'))
    evals, Q = eigen(S)
    X = X*Q
    if verbose_last || verbose
        @printf("\nEigenvalues from diag_LOBPCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%d] = %18.10f, devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end
    return evals, X
end
