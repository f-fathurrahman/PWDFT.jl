"""
Locally-optimal block preconditioned conjugate gradient method for
finding `Hamiltonian` eigenstates and eigenvalues.

Based on code by Knyazev.

TODO: Add references
"""
function diag_lobpcg( Ham::Hamiltonian, X0::Array{ComplexF64,2};
                      tol=1e-6, NiterMax=200, verbose=false,
                      verbose_last=false, Nstates_conv=0 )

    pw = Ham.pw
    # get size info
    Nbasis = size(X0)[1]
    Nstates = size(X0)[2]

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    # orthonormalize the initial wave functions.
    #X = ortho_gram_schmidt(X0)  # normalize (again) XXX ?
    X = copy(X0)

    HX = op_H( Ham, X )

    nconv = 0
    iter = 1
    resnrm = ones(Nstates,1)

    P = zeros(ComplexF64,Nbasis,Nstates)
    HP = zeros(ComplexF64,Nbasis,Nstates)

    sum_evals = 0.0
    sum_evals_old = 0.0
    conv = 0.0

    IS_CONV = false
    
    for iter = 1:NiterMax

        # Rayleigh quotient
        S = X'*HX
        lambda = real(eigvals(S))  #
        #
        # Check for convergence
        #
        sum_evals = sum(lambda)
        conv = abs(sum_evals - sum_evals_old)/Nstates
        sum_evals_old = sum_evals
        R = HX - X*S

        nconv = length( findall( resnrm .< tol ) )
        if nconv >= Nstates_conv
            IS_CONV = true
            if verbose
                println("LOBPCG convergence: nconv in iter ", iter)
            end
            break
        end

        if verbose
            @printf("LOBPCG iter = %8d, nconv=%d, %18.10e\n", iter, nconv, conv)
        end
        

        for j = 1:Nstates
            resnrm[j] = norm( R[:,j] )
        end
        #
        # apply preconditioner
        W = Kprec( Ham.ik, pw, R )
        #
        # nlock == 0
        #
        HW = op_H( Ham, W )
        #
        C  = W'*W
        #C = 0.5*( C + C' )
        C = Hermitian(C)
        #
        R  = (cholesky(C)).U
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
            C = Hermitian(C)
            R = (cholesky(C)).U
            P = P/R
            HP = HP/R
        else
            P  = copy(W)
            HP = copy(HW)
        end

        iter = iter + 1
    end

    if !IS_CONV
        @printf("\nWARNING: LOBPCG is not converged after %d iterations\n", NiterMax)
    end

    S = X'*HX
    S = (S+S')/2
    lambda, Q = eigen(S)
    X = X*Q
    if verbose_last
        for j = 1:Nstates
            @printf("eigval[%2d] = %18.10f, resnrm = %18.10e\n", j, lambda[j], resnrm[j] )
        end
    end
    #@printf("LOBPCG converge: %d iter\n", iter)
    return lambda, X
end
