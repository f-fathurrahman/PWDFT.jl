function linmin_armijo!( Ham,
    psiks_orig, g, d, E_orig; α0 = 1.0, τ=0.1, c=0.1 )

    psiks = deepcopy(psiks_orig)
    m = dot_BlochWavefunc( g, d )
    @printf("m = %18.10f\n", m)

    t = -c*m
    α = α0

    for itry = 1:3
        do_step!( psiks, α, d )
        E_trial = calc_energies_only!( Ham, psiks )
        dE = E_orig - E_trial
        @printf("α = %e, E_trial = %e, dE = %e, α*t = %e\n", α, E_trial, dE, α*t)
        #if dE < 0.0
        #    @printf("ERROR: α is not reducing E_trial\n")
        #end
        if dE >= α*t
            @printf("Condition is satisfied at α = %e\n", α)
            return true, α
        end
        α = τ*α
        psiks = deepcopy(psiks_orig)
    end
    @printf("**WARNING: Condition is not satisfied, returning: α=%e\n", α)
    return true, α

end