function linmin_debug!( Ham, psiks;
    etot_conv_thr=1e-7, skip_initial_diag=false, startingrhoe=:gaussian
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

    Etot = calc_energies_grad!( Ham, psiks, g, Kg )
    #println("Starting Etot = ", Etot)
    #println("dot_BlochWavefunc(g,g) = ", dot_BlochWavefunc(g,g))

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, psiks )

    N_α_adjust_max = 3
    αt_start = 1.0 # This should be a parameter
    αt = αt_start
    α = αt
    αt_min = 1e-8 #1e-10
    αt_reduceFactor = 0.1 # should less than 1
    αt_increaseFactor = 3.0
    updateTestStepSize = true

    β = 0.0
    gKnorm = 0.0
    gKnormPrev = 0.0
    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    gKnorm = dot_BlochWavefunc(g, Kg)

    β = 0.0

    # Update search direction
    for i in 1:Nkspin
        d[i] = -Kg[i] + β*d[i]
    end

    constrain_search_dir!( d, psiks )

    αPrev = 0.0

    psiks_t = zeros_BlochWavefunc(Ham)

    println()

    test_undoing!(Ham, psiks, d)

    #αt = 1.0
    ## Try the test step
    #for i in 1:Nkspin
    #    psiks_t[i] = psiks[i] + αt*d[i]
    #end
    #E_trial = calc_energies_only!( Ham, psiks_t )
    ##E_trial = calc_energies_grad!( Ham, psiks, g, Kg )
    #@printf("αt = %18.10e, E_trial = %18.10f\n", αt, E_trial)

    #αt = 2.0
    ## Try the test step
    #for i in 1:Nkspin
    #    psiks_t[i] = psiks[i] + αt*d[i]
    #end
    #E_trial = calc_energies_only!( Ham, psiks_t )
    ##E_trial = calc_energies_grad!( Ham, psiks, g, Kg )
    #@printf("αt = %18.10e, E_trial = %18.10f\n", αt, E_trial)

end


function test_undoing!( Ham, psiks, d )

    Nkspin = length(psiks)
    
    E0 = calc_energies_only!( Ham, psiks )
    println("E0 = ", E0)

    αt = 1.0
    # Try the test step
    for i in 1:Nkspin
        psiks[i] = psiks[i] + αt*d[i]
    end
    #E_trial = calc_energies_only!( Ham, psiks )
    E_trial = calc_energies_no_modify!( Ham, psiks )
    @printf("E_trial = %18.10f\n", E_trial)

    for i in 1:Nkspin
        psiks[i] = psiks[i] - αt*d[i]
    end
    #E_trial = calc_energies_only!( Ham, psiks )
    E_trial = calc_energies_no_modify!( Ham, psiks )
    @printf("Undoing = %18.10f\n", E_trial)
end
