"""
Based on code by Knyazev.
"""
function diag_lobpcg( Ham::PWHamiltonian, X0::Array{Complex128,2};
                      tol=1e-5, tol_avg=1e-7, maxit=200, verbose=false,
                      verbose_last=false )

    pw = Ham.pw
    # get size info
    Nstates = size(X0)[2]

    # orthonormalize the initial wave functions.
    X = ortho_gram_schmidt(X0)  # normalize (again)

    HX = op_H( Ham, X )

    nconv = 0
    iter = 1
    resnrm = ones(Nstates,1)

    sum_evals = 0.0
    sum_evals_old = 0.0
    conv = 0.0
    while iter <= maxit && nconv < Nstates
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
        if verbose
            @printf("LOBPCG iter = %8d, %18.10e\n", iter, conv)
        end
        if conv <= tol_avg
            if verbose
                @printf("LOBPCG convergence: tol_avg\n")
            end
            break
        end
        #
        for j = 1:Nstates
            resnrm[j] = norm( R[:,j] )
        end
        #
        # apply preconditioner
        W = Kprec(pw,R)
        #
        # nlock == 0
        #
        HW = op_H( Ham, W )
        #
        C  = W'*W
        C = ( C + C' )/2
        R  = chol(C)
        W  = W/R
        HW = HW/R
        #
        Q  = [X W]
        HQ = [HX HW]
        if iter > 1
            Q  = [Q P]
            HQ = [HQ HP]
        end

        T = Q'*(HQ); T = (T+T')/2;
        G = Q'*Q; G = (G+G')/2;

        sd, S = eig( T, G ) # evals, evecs
        U = S[:,1:Nstates]
        X = Q*U
        HX = HQ*U
        if iter > 1
            set2 = Nstates  +1:2*Nstates
            set3 = 2*Nstates + 1:3*Nstates
            P  = W*U[set2,:]  + P*U[set3,:]
            HP = HW*U[set2,:] + HP*U[set3,:]
            C = P'*P
            C = (C + C')/2
            R = chol(C)
            P = P/R
            HP = HP/R
        else
            P  = copy(W)
            HP = copy(HW)
        end

        iter = iter + 1
    end

    S = X'*HX
    S = (S+S')/2
    lambda, Q = eig(S)
    X = X*Q
    if verbose_last
        for j = 1:Nstates
            @printf("eigval[%2d] = %18.10f, resnrm = %18.10e\n", j, lambda[j], resnrm[j] )
        end
    end
    #@printf("LOBPCG converge: %d iter\n", iter)
    return lambda, X
end
