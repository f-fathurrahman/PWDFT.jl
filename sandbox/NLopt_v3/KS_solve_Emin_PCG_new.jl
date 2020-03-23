function KS_solve_Emin_PCG_new!( Ham, psiks;
    etot_conv_thr=1e-6, skip_initial_diag=false, startingrhoe=:gaussian
)

    Nkspin = length(psiks)

    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    gPrev = zeros_BlochWavefunc(Ham)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates
    Rhoe = zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end

    update!(Ham, Rhoe)

    evals = zeros(Nstates,Nkspin)
    if !skip_initial_diag
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )
    end

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )
    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    Etot = calc_energies_grad!( Ham, psiks, g, Kg, skip_ortho=true )
    println("Starting Etot = ", Etot)

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, psiks )

    gPrevUsed = true

    αt_start = 1.0 # This should be a parameter
    αt = αt_start
    α = αt
    αt_min = 1e-10
    αt_reduceFactor = 0.1 # should less than 1
    αt_increaseFactor = 3.0
    updateTestStepSize = true

    β = 0.0
    gKnorm = 0.0
    gKnormPrev = 0.0
    force_grad_dir = true

    NiterMax = 50

    Etot_old = Etot
    Nconverges = 0
    
    for iter in 1:NiterMax
        gKnorm = dot_BlochWavefunc(g, Kg)

        β = 0.0
        
        if !force_grad_dir
            
            dotgd = dot_BlochWavefunc(g, d)
            if gPrevUsed
                dotgPrevKg = dot_BlochWavefunc(gPrev, Kg)
            else
                dotgPrevKg = 0.0
            end

            β = (gKnorm - dotgPrevKg)/gKnormPrev
        end

        println("\nBegin iter = ", iter)
        println("β raw = ", β)
        if β < 0.0
            println("Resetting β")
            β = 0.0
        end

        force_grad_dir = false

        # Check convergence here ....

        # No convergence yet, continuing ...
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = gKnorm

        # Update search direction
        d = -Kg + β*d

        constrain_search_dir!( d, psiks )

        # Line minimization
        #linmin_success, Etot, α, αt = linmin_quad!( Ham, psiks, g, Kg, d, α, αt, Etot )
        #println("linmin_success = ", linmin_success)
        #@printf("α = %f, αt = %f\n", α, αt)
1
        linmin_success = true
        Etot = linmin_grad!( Ham, psiks, g, Kg, d )

        if linmin_success
            if updateTestStepSize
                αt = α
                if αt < αt_min    # bad step size
                    @printf("Bad step size: αt is reset to αt_start: %f\n", αt_start)
                    αt = αt_start # make sure next test step size is not too bad
                end
            end
        else
            # linmin failed:
            @printf("Undoing step.\n")
            #step(d, -alpha)
            for i in 1:Nkspin
                psiks[i] = psiks[i] - α*d[i]
            end

            Etot = calc_energies_grad!( Ham, psiks, g, Kg )
            println("Etot after undoing step: ", Etot)
            
            if β > 0.0
                # Failed, but not along the gradient direction:
                @printf("Step failed: resetting search direction.\n")
                forceGradDirection = true # reset search direction
            else
                # Failed along the gradient direction
                @printf("Step failed along negative gradient direction.\n")
                @printf("Probably at roundoff error limit. (Should stop)\n")
                #return
            end
        end
        
        diffE = abs(Etot - Etot_old)
        @printf("Emin_PCG_new step %8d = %18.10f %10.7e\n", iter, Etot, diffE)
        
        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end
        
        if Nconverges >= 2
            @printf("\nEmin_PCG is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot

    end


    println("Leaving KS_solve_Emin_PCG_new")

end


function constrain_search_dir!( d, psiks )
    Nkspin = length(psiks)
    for i in 1:Nkspin
        d[i] = d[i] - psiks[i] * ( psiks[i]' * d[i] )
    end
    return
end

function dot_BlochWavefunc(x::BlochWavefunc, y::BlochWavefunc)
    Nkspin = length(x)    
    res = 2.0
    for i in 1:Nkspin
        res = res + real( dot(x[i], y[i]) )*2.0
    end
    return res
end

