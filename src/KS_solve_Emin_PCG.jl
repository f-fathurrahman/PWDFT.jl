#
# Ham.potentials.V_Ps_loc should be initialized
# Ham.energies.NN should be calculated if needed
#
function KS_solve_Emin_PCG!( Ham::PWHamiltonian;
                             α_t = 3e-5, NiterMax=1000, verbose=false,
                             I_CG_BETA = 2 )

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Ns = pw.Ns
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    #
    # Random guess of wave function
    #
    srand(1234)
    psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    rhoe = calc_rhoe( pw, Focc, psi )
    update!(Ham, rhoe)

    #
    # Variabls for PCG
    #
    d = zeros(Complex128, Ngwx, Nstates)
    g_old = zeros(Complex128, Ngwx, Nstates)
    d_old = zeros(Complex128, Ngwx, Nstates)
    Kg = zeros(Complex128, Ngwx, Nstates)
    Kg_old = zeros(Complex128, Ngwx, Nstates)
    β        = 0.0
    Etot_old = 0.0

    # Calculate energy at this psi
    Energies = calc_energies(Ham, psi)
    Ham.energies = Energies

    Etot     = Energies.Total

    @printf("Initial Etot = %18.10f\n", Etot)


    for iter = 1:NiterMax

        g = calc_grad( Ham, psi)
        Kg = Kprec(pw,g)

        #println("sum g  = ", sum(g))
        #println("sum Kg = ", sum(Kg))

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
        rhoe = calc_rhoe( pw, Focc, psi )

        update!(Ham, rhoe)

        Energies = calc_energies( Ham, psi )
        Ham.energies = Energies
        Etot = Energies.Total
        diff = abs(Etot-Etot_old)

        @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        if diff < 1e-6
            @printf("CONVERGENCE ACHIEVED\n")
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot
    end

    psi = ortho_gram_schmidt(psi)
    Hr = psi' * op_H( Ham, psi )
    evals, evecs = eig(Hr)
    psi = psi*evecs

    return real(evals), psi
end
