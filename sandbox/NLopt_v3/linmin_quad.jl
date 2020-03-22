# return succeed or not
# alpha return E
# modify psiks
function linmin_quad!( Ham, psiks, g, Kg, d, αt, α, E )
    αPrev = 0.0
    Eorig = E
    gdotd = dot_BlochWavefunc(g,d) # directional derivative at starting point

    # These should be parameters
    N_α_adjust_max = 3
    αt_min = 1e-10
    αt_reduceFactor = 0.1 # should less than 1
    αt_increaseFactor = 3.0

    Nkspin = length(psiks)

    for s in 1:N_α_adjust_max
        
        if αt < αt_min
            println("αt below threshold. Quitting step.")
            alpha = αPrev
            return false, E
        end

        # Try the test step
        #step(d, αt - αPrev)
        # axpy(alpha, rotExists ? dir.C[q]*rotPrevC[q] : dir.C[q], eVars.C[q]);
        for i in 1:Nkspin
            psiks[i] = psiks[i] + (αt - αPrev)*d[i]
        end
        αPrev = αt
        
        ET = calc_energies_only!( Ham, psiks )
        
        # Check if step crossed domain of validity of parameter space:
        #if(!std::isfinite(ET)) SKIPPED

        # Predict step size:
        α = 0.5 * αt^2 *gdotd / (αt * gdotd + E - ET)

        # Check reasonableness of predicted step size:
        
        if α < 0    
            # Curvature has the wrong sign
            # That implies ET < E, so accept step for now, and try descending further next time
            αt = αt * αt_increaseFactor;
            @printf("Wrong curvature in test step, increasing αt to %le.\n", αt)
            
            E = calc_energies_grad!( Ham, psiks, g, Kg )
            return true, E  # also return new energy
        end
        
        if α/αt > αt_increaseFactor
            αt = αt * αt_increaseFactor
            @printf("Predicted α/αt > %lf, increasing alphaT to %le.\n", αt_increaseFactor, αt)
            continue
        end

        if αt/α < αt_reduceFactor
            αt = αt * αt_reduceFactor;
            @printf("Predicted α/αt < %lf, reducing alphaT to %le.\n", αt_reduceFactor, αt)
            continue
        end
        
        # Successful test step:
        break
    end


    if !isfinite(E)
        @printf("Test step failed %d times. Quitting step.\n", N_α_adjust_max)
        α = αPrev
        return false, E # E should not be used
    end


    # Actual step:
    for s in 1:N_α_adjust_max
        
        # Try the step:
        #obj.step(d, alpha-αPrev); αPrev=alpha;
        for i in 1:Nkspin
            psiks[i] = psiks[i] + (α - αPrev)*d[i]
        end
        αPrev = α
        
        E = calc_energies_grad!( Ham, psiks, g, Kg )
        
        if !isfinite(E)
            α = α * αt_reduceFactor;
            @printf("Step failed: E = %le, reducing alpha to %le.\n", E, α)
            continue
        end
        
        if E > Eorig
            α = α * αt_reduceFactor
            @printf("Step increased by: %le, reducing alpha to %le.\n", E-Eorig, α)
            continue
        end
        
        # Step successful:
        break
    end
    
    if !isfinite(E) || (E > Eorig)
        @printf("Step failed to reduce after %d attempts. Quitting step.\n", N_α_adjust_max)
        return false, E
    end
    
    return true, E

end