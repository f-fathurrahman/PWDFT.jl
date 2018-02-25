function KS_solve_SCF!( Ham::PWHamiltonian, Nstates::Int64;
                       β = 0.5, NiterMax=100, verbose=false,
                       check_rhoe_after_mix=false )

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints
    Focc = Ham.focc

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


    Etot_old = 0.0
    rhoe_old = copy(rhoe)

    λ = zeros(Float64,Nstates)

    for iter = 1:NiterMax

        λ, psi = diag_lobpcg( Ham, psi, verbose_last=false )

        #
        rhoe_new = calc_rhoe( pw, Focc, psi )
        diffRho = norm(rhoe_new - rhoe)/Npoints
        #
        rhoe = β*rhoe_new[:] + (1-β)*rhoe[:]

        if check_rhoe_after_mix
            integRhoe = sum(rhoe)*ΔV
            @printf("After mixing: integRho = %18.10f\n", integRhoe)
        end

        update!( Ham, rhoe )

        # Calculate energies
        Energies = calc_energies( Ham, psi )
        Ham.energies = Energies
        Etot = Energies.Total
        diffE = abs( Etot - Etot_old )

        #
        @printf("SCF: %8d %18.10f %18.10e %18.10e\n", iter, Etot, diffE, diffRho )

        if diffE < 1e-6
            @printf("SCF is is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end
        #
        Etot_old = Etot
    end

    return λ, psi

end
