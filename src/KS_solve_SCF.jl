function KS_solve_SCF!( Ham::PWHamiltonian, Nstates::Int64;
                       β = 0.5, NiterMax=100, verbose=false,
                       check_rhoe_after_mix=false,
                       update_psi="LOBPCG",
                       cheby_degree=8 )

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

    const ETHR_EVALS_LAST = 1e-6

    ethr = 0.1

    #
    # For mixing
    #
    MIXDIM = 4
    df = zeros(Float64,Npoints,MIXDIM)
    dv = zeros(Float64,Npoints,MIXDIM)


    for iter = 1:NiterMax

        if update_psi == "LOBPCG"
            λ, psi = diag_lobpcg( Ham, psi, verbose_last=false )

        elseif update_psi == "PCG"

            if iter == 1
                ethr = 0.1
            elseif iter == 2
                ethr = 0.01
            else
                ethr = ethr/5.0
                ethr = max( ethr, ETHR_EVALS_LAST )
            end

            λ, psi = diag_Emin_PCG( Ham, psi, TOL_EBANDS=ethr )

        elseif update_psi == "CheFSI"

            ub, lb = get_ub_lb_lanczos( Ham, Nstates*2 )

            psi = chebyfilt( Ham, psi, cheby_degree, lb, ub)
            psi = ortho_gram_schmidt(psi)

        end

        #
        rhoe_new = calc_rhoe( pw, Focc, psi )
        diffRho = norm(rhoe_new - rhoe)/Npoints

        #rhoe = β*rhoe_new[:] + (1-β)*rhoe[:]

        rhoe = andersonmix!( rhoe, rhoe_new, β, df, dv, iter, MIXDIM )

        for ip = 1:Npoints
            if rhoe[ip] < 1e-12
                rhoe[ip] = 1e-12
            end
        end

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

    if update_psi == "CheFSI"
        Hr = psi'*op_H( Ham, psi )
        λ = real(eigvals(Hr))
    end

    return λ, psi

end
