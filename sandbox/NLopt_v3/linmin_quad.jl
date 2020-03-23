# return succeed or not
# modify Ham, psiks, g, Kg
# return success_or_not, E, α, αt
function linmin_quad!( Ham, psiks, g, Kg, d, α, αt, E )

    Nkspin = length(psiks)

    αPrev = 0.0
    E_orig = E
    gdotd = dot_BlochWavefunc(g,d) # directional derivative at starting point
    #for i in 1:Nkspin
    #    println("dot g d: ", dot(g[i], d[i]))
    #end
    println("gdotd = ", gdotd)

    if gdotd >= 0.0
        @printf("Bad step direction: g.d = %f > 0.0\n", gdotd)
        α = αPrev
        return false, E, α, αt
    end

    # These should be parameters
    N_α_adjust_max = 3
    αt_min = 1e-8 #1e-10
    αt_reduceFactor = 0.1 # should less than 1
    αt_increaseFactor = 3.0

    E_trial = 0.0
    
    #psiks_t = deepcopy(psiks)

    for s in 1:N_α_adjust_max
        
        if αt < αt_min
            println("αt below threshold. Quitting step.")
            α = αPrev
            return false, E, α, αt
        end

        # Try the test step
        #step(d, αt - αPrev)
        # axpy(alpha, rotExists ? dir.C[q]*rotPrevC[q] : dir.C[q], eVars.C[q]);
        for i in 1:Nkspin
            psiks[i] = psiks[i] + (αt - αPrev)*d[i]
        end
        αPrev = αt
        #E_trial = calc_energies_only!( Ham, psiks_t )
        E_trial = calc_energies_grad!( Ham, psiks, g, Kg )
        println("E_trial = ", E_trial)
        
        # Check if step crossed domain of validity of parameter space:
        if !isfinite(E_trial)
            αt = αt * αt_reduceFactor
            @printf("Test step failed, E_trial = %le, reducing αt to %le.\n", E_trial, αt)
            continue
        end


        # Predict step size:
        α = 0.5 * αt^2 *gdotd / (αt * gdotd + E - E_trial)

        # Check reasonableness of predicted step size:
        
        if α < 0
            # Curvature has the wrong sign
            # That implies ET < E, so accept step for now, and try descending further next time
            αt = αt * αt_increaseFactor;
            @printf("Wrong curvature in test step, increasing αt to %le.\n", αt)
            
            E = calc_energies_grad!( Ham, psiks, g, Kg )
            @printf("Before return: α = %f, αt = %f\n", α, αt)
            @printf("linmin: dEtot = %18.10e\n", E_orig - E)
            return true, E, α, αt  # also return new energy
        end
        
        if α/αt > αt_increaseFactor
            αt = αt * αt_increaseFactor
            @printf("Predicted α/αt > %lf, increasing αt to %le.\n", αt_increaseFactor, αt)
            continue
        end

        if αt/α < αt_reduceFactor
            αt = αt * αt_reduceFactor
            @printf("Predicted α/αt < %lf, reducing αt to %le.\n", αt_reduceFactor, αt)
            continue
        end
        
        # Successful test step:
        break
    end


    if !isfinite(E)
        @printf("Test step failed %d times. Quitting step.\n", N_α_adjust_max)
        α = αPrev
        return false, E, α, αt # E should not be used
    end


    # Actual step:
    α_total = 0.0
    for s in 1:N_α_adjust_max
        
        # Try the step:
        #obj.step(d, alpha-αPrev); αPrev=alpha;
        for i in 1:Nkspin
            psiks[i] = psiks[i] + (α - αPrev)*d[i]
        end
        α_total = α_total + (α - αPrev)
        αPrev = α
        
        E = calc_energies_grad!( Ham, psiks, g, Kg )
        @printf("linmin actual step: α = %18.10e E = %18.10f\n", α, E)
        
        if !isfinite(E)
            α = α * αt_reduceFactor;
            @printf("Step failed: E = %le, reducing α to %le.\n", E, α)
            continue
        end
        
        if E > E_orig
            α = α * αt_reduceFactor
            @printf("Step increased by: %le, reducing α to %le.\n", E - E_orig, α)
            continue
        end
        
        # Step successful:
        break
    end
    
    if !isfinite(E) || (E > E_orig)
        @printf("Step failed to reduce after %d attempts. Quitting step.\n", N_α_adjust_max)
        #println("α_total = ", α_total)
        return false, E, α, αt
    end
    
    dE_trial = E_trial - E
    if dE_trial < 0.0
        @printf("*** WARNING linmin: E_trial - E is negative: %18.10e\n", dE_trial)
    end
    @printf("linmin: dEtot = %18.10e\n", E_orig - E)

    return true, E, α, αt

end