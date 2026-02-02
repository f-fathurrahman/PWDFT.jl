function diag_Emin_PCG( Ham::HamiltonianGamma, X0::Array{ComplexF64,2};
                        tol=1e-5, NiterMax=100, verbose=false,
                        verbose_last=false, Nstates_conv=0,
                        tol_ebands=1e-4,
                        α_t=3e-5, I_CG_BETA = 2 )

    X = copy(X0)
    evals = diag_Emin_PCG!( Ham, X, tol=tol, NiterMax=NiterMax, verbose=verbose,
                            verbose_last=verbose_last, Nstates_conv=Nstates_conv,
                            tol_ebands=tol_ebands,
                            α_t=α_t)                          
    return evals, X
end


function diag_Emin_PCG!( Ham::HamiltonianGamma, psis::BlochWavefuncGamma;
                         tol=1e-5, NiterMax=100, verbose=false,
                         verbose_last=false, Nstates_conv=0,
                         tol_ebands=1e-4,
                         α_t=3e-5)
    
    pw = Ham.pw
    electrons = Ham.electrons
    Nspin = electrons.Nspin_wf
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates
    evals = zeros(Float64,Nstates,Nspin)
    for ispin = 1:Nspin
        Ham.ispin = ispin
        evals[:,ispin] =
        diag_Emin_PCG!( Ham, psis.data[ispin], tol=tol, NiterMax=NiterMax, verbose=verbose,
                       verbose_last=verbose_last, Nstates_conv=Nstates_conv,
                       tol_ebands=tol_ebands,
                       α_t=α_t )
    end
    return evals
end


function diag_Emin_PCG!( Ham::HamiltonianGamma, X::Array{ComplexF64,2};
                         tol=1e-5, NiterMax=100, verbose=false,
                         verbose_last=false, Nstates_conv=0,
                         tol_ebands=1e-4,
                         α_t=3e-5 )

    pw = Ham.pw
    Focc = Ham.electrons.Focc

    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(X)[2]
    Ngw = pw.gvecw.Ngw

    if Nstates_conv == 0
        Nstates_conv = Nstates
    end

    #
    # Variabls for PCG
    #
    d = zeros(ComplexF64, Ngw, Nstates)
    g = zeros(ComplexF64, Ngw, Nstates)
    gPrev = zeros(ComplexF64, Ngw, Nstates)
    d_old = zeros(ComplexF64, Ngw, Nstates)
    Kg = zeros(ComplexF64, Ngw, Nstates)
    Xc = zeros(ComplexF64, Ngw, Nstates)
    gt = zeros(ComplexF64, Ngw, Nstates)
    Hr = zeros(ComplexF64, Nstates, Nstates)
    
    β = 0.0
    
    Ebands_old = 0.0

    calc_grad_evals!(Ham, X, g, Hr)
    
    evals = eigvals(Hermitian(Hr))
    Ebands = sum(evals)
    
    evals_old = copy(evals)
    devals = ones(Nstates)

    IS_CONVERGED = false
    gKnormPrev = 1.0

    for iter = 1:NiterMax

        Kprec!(pw, g, Kg)

        gKnorm = 2*real(dot_gamma(g, Kg))
        if iter != 1
            dotgd = 2*real(dot_gamma(g, d))
            dotgPrevKg = 2*real(dot_gamma(gPrev, Kg))
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        end
        if β < 0.0
            β = 0.0
        end

        # Save old gradient
        gPrev[:,:] = g[:,:]
        gKnormPrev = gKnorm

        # Update search direction
        d = -Kg + β * d_old

        Xc[:,:] = X + α_t*d
        ortho_sqrt_gamma!(Xc)

        # Note that Ham is not updated

        calc_grad_evals!(Ham, Xc, gt, Hr)

        denum = 2.0*real(dot_gamma(g-gt, d))
        if denum != 0.0
            α = abs(α_t * 2.0*real(dot_gamma(g, d))/denum)
        else
            α = 0.0
        end

        # Update wavefunction
        X[:,:] = X + α*d
        ortho_sqrt_gamma!(X)

        calc_grad_evals!(Ham, X, g, Hr)
        
        evals = eigvals(Hermitian(Hr))
        Ebands = sum(evals)

        devals = abs.( evals - evals_old )
        evals_old = copy(evals)
        
        nconv = length( findall( devals .< tol ) )

        diffE = abs(Ebands-Ebands_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diffE)
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
        if diffE <= tol_ebands*Nstates
            IS_CONVERGED = true
            if verbose || verbose_last
                @printf("Convergence is achieved based on tol_ebands*Nstates\n")
            end
            break
        end

        d_old = copy(d)
        Ebands_old = Ebands
    end

    if !IS_CONVERGED
        @printf("\nWARNING: diag_Emin_PCG is not converged after %d iterations\n", NiterMax)
    end

    ortho_sqrt_gamma!(X)
    
    Hr = overlap_gamma(X, op_H(Ham, X))
    evals, evecs = eigen(Hermitian(Hr))
    X[:,:] = X*evecs

    if verbose_last || verbose
        @printf("\nEigenvalues from diag_Emin_PCG:\n\n")
        for ist = 1:Nstates
            @printf("evals[%3d] = %18.10f devals = %18.10e\n", ist, evals[ist], devals[ist] )
        end
    end

    return evals
end
