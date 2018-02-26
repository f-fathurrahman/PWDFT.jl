#
# Ham.potentials.V_Ps_loc should be initialized
# Ham.energies.NN should be calculated if needed
#
function diag_Emin_PCG( Ham::PWHamiltonian, X0::Array{Complex128,2};
                             α_t = 3e-5, NiterMax=200, verbose=false,
                             I_CG_BETA = 1, TOL_EBANDS=1e-7 )


    pw = Ham.pw
    Focc = Ham.focc

    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(X0)[2]
    Ngwx = pw.gvecw.Ngwx

    psi = ortho_gram_schmidt(X0)

    #
    # Variabls for PCG
    #
    d = zeros(Complex128, Ngwx, Nstates)
    g_old = zeros(Complex128, Ngwx, Nstates)
    d_old = zeros(Complex128, Ngwx, Nstates)
    Kg = zeros(Complex128, Ngwx, Nstates)
    Kg_old = zeros(Complex128, Ngwx, Nstates)
    β        = 0.0
    Ebands_old = 0.0

    Hr = psi' * op_H( Ham, psi )
    Ebands = real(sum(eigvals(Hr)))

    @printf("\nInitial Ebands = %18.10f\n", Ebands)

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
            @printf("β is smaller than 0, setting it to zero\n")
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

        Hr = psi' * op_H( Ham, psi )
        Ebands = real(sum(eigvals(Hr)))

        diff = abs(Ebands-Ebands_old)

        @printf("CG step %8d = %18.10f %10.7e\n", iter, Ebands, diff)
        if diff < TOL_EBANDS
            @printf("CONVERGENCE ACHIEVED\n")
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Ebands_old = Ebands
    end

    psi = ortho_gram_schmidt(psi)
    Hr = psi' * op_H( Ham, psi )
    evals, evecs = eig(Hr)
    psi = psi*evecs

    return real(evals), psi
end
