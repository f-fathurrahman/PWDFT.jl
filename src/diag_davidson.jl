function diag_davidson( Ham::PWHamiltonian, X0::Array{ComplexF64,2};
                        tol=1e-5, tol_avg=1e-7, NiterMax=200, verbose=false,
                        verbose_last=false )
      

    # get size info
    Nstates = size(X0)[2]
    Nbasis  = size(X0)[1]

    if Nstates <= 0
        @printf("diag_davidson requires at least one initial wave function!\n");
        return
    end
    
    # orthonormalize the initial wave functions.
    X = ortho_gram_schmidt(X0)  # normalize (again)?

    evals    = zeros(Float64, Nstates)
    R        = zeros(ComplexF64, Nbasis, Nstates)
    Hred     = zeros(ComplexF64, 2*Nstates, 2*Nstates)
    Sred     = zeros(ComplexF64, 2*Nstates, 2*Nstates)
    res      = zeros(Float64, Nstates)
    res_norm = zeros(Float64, Nstates)
    devals   = zeros(Float64, Nstates)

    HX = op_H( Ham, X )

    # Initial eigenvalues
    for ic = 1:Nstates
        for ig = 1:Nbasis
            #evals[ic] = evals[ic] + real( conj(X[ig,ic]) * HX[ig,ic] )
            evals[ic] = real( conj(X[ig,ic]) * HX[ig,ic] )
        end
    end

    # Calculate residuals
    for ic = 1:Nstates
        for ig = 1:Nbasis
            R[ig,ic] = evals[ic]*X[ig,ic] - HX[ig,ic]
        end
    end

    for ic = 1:Nstates
        res[ic] = 0.0
        for ig = 1:Nbasis
            res[ic] = res[ic] + real( R[ig,ic] * conj(R[ig,ic]) )
        end
        res[ic] = sqrt( res[ic] )
    end

    EPS = eps()

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates

    sum_evals = sum(evals)
    sum_evals_old = sum_evals

    for iter = 1:NiterMax

        res_norm[:] .= 1.0

        for ic = 1:Nstates
            if EPS < res[ic]
                res_norm[ic] = 1.0/res[ic]
            end
        end

        for ic = 1:Nstates
            for ig = 1:Nbasis
                R[ig,ic] = res_norm[ic] * R[ig,ic]
            end
        end

        R = Kprec(Ham.ik, Ham.pw, R)

        HR = op_H(Ham, R)

        # FIXME: Pull this outside the loop?
        if iter == 1
            Hred[set1,set1] = X' * HX
        else
            Hred[1:Nstates,1:Nstates] .= 0.0 + im*0.0
            for ic = 1:Nstates
                Hred[ic,ic] = evals[ic] + im*0.0 # use diagm ?
            end
        end

        Hred[set1,set2] = X' * HR
        Hred[set2,set2] = R' * HR
        Hred[set2,set1] = Hred[set1,set2]'

        Sred[set1,set1] = Matrix(Diagonal(ones(ComplexF64,Nstates)))
        Sred[set1,set2] = X' * R
        Sred[set2,set2] = R' * R
        Sred[set2,set1] = Sred[set1,set2]'

        Hred = (Hred + Hred')/2.0
        Sred = (Sred + Sred')/2.0

        λ_red, X_red = eigen( Hermitian(Hred), Hermitian(Sred) )

        evals = λ_red[1:Nstates]

        #for ic = 1:2*Nstates
        #  print(ic); print(" "); println(λ_red[ic])
        #end

        X  = X  * X_red[set1,set1] + R  * X_red[set2,set1]
        HX = HX * X_red[set1,set1] + HR * X_red[set2,set1]

        # Calculate residuals
        for ic = 1:Nstates
            for ig = 1:Nbasis
                R[ig,ic] = evals[ic]*X[ig,ic] - HX[ig,ic]
            end
        end

        for ic = 1:Nstates
            res[ic] = 0.0
            for ig = 1:Nbasis
                res[ic] = res[ic] + real( R[ig,ic] * conj(R[ig,ic]) )
            end
            res[ic] = sqrt( res[ic] )
        end
        
        sum_evals = sum(evals)
        diffSum = abs(sum_evals - sum_evals_old)
        if verbose
            @printf("iter %d tol_avg = %18.10e\n", iter, abs(diffSum))
        end
        if diffSum < tol_avg
            if verbose || verbose_last
                @printf("Davidson convergence: tol_avg in iter = %d\n", iter)
            end
            break
        end
        sum_evals_old = sum_evals
    end

    if verbose_last
        for j = 1:Nstates
            @printf("eigval[%2d] = %18.10f, resnrm = %18.10e\n", j, evals[j], res[j] )
        end
    end


    return evals, X

end
