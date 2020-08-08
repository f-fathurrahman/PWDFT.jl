using NLsolve

"""
NLsolve driver.
Adapted from DFTK.jl
"""
function scf_NLsolve_solver( m=5; kwargs... )

    function fix_point_solver( f, x0, tol, NiterMax )
        res = NLsolve.nlsolve(x -> f(x) - x, x0;
            m=m, method=:anderson, xtol=tol, ftol=0.0, show_trace=true,
            iterations=NiterMax, kwargs... )
        ( fixpoint=res.zero, converged=NLsolve.converged(res) )
    end

    return fix_point_solver

end


# A draft implementation of using NLsolve for
function KS_solve_SCF_NLsolve!( Ham::Hamiltonian; use_smearing=false, kT=1e-3 )

    Random.seed!(1234)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Nstates_occ = Ham.electrons.Nstates_occ

    Rhoe = zeros(Float64, Npoints, Nspin)
    psiks = rand_BlochWavefunc(Ham)
    
    #calc_rhoe!( Ham, psiks, Rhoe )
    Rhoe[:,1] = guess_rhoe( Ham )

    # calculate E_NN and PspCore energies
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    evals = zeros( Float64, Nstates, Nkspin )

    ethr_last = 1e-5
    betamix = 1.0 # use Rhoe_out directly

    Entropy = 0.0
    E_fermi = 0.0
    Nelectrons = Ham.electrons.Nelectrons
    Focc = copy(Ham.electrons.Focc) # make sure to use the copy
    wk = Ham.pw.gvecw.kpoints.wk

    function density_map( Rhoe_in )
        update!( Ham, Rhoe_in )
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr_last,
                      Nstates_conv=Nstates_occ )
        
        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        Rhoe_out = calc_rhoe( Ham, psiks )
        # simple mixing here?
        #Rhoe_next = betamix*Rhoe_out + (1 - betamix)*Rhoe_in
        return Rhoe_out
    end

    NiterMax = 100

    println("\nNLsolve starting: \n")

    fpres = NLsolve.nlsolve(x -> density_map(x) - x, Rhoe,
            m=5, method=:anderson, xtol=1e-6, ftol=0.0, show_trace=true,
            iterations=NiterMax )
    
    println("\nNLsolve is converged in iter: ", fpres.iterations)
    
    #println(fieldnames(typeof(fpres)))
    # fieldnames: (:method, :initial_x, :zero, :residual_norm, :iterations, :x_converged,
    # :xtol, :f_converged, :ftol, :trace, :f_calls, :g_calls)

    update!( Ham, fpres.zero )

    println()
    println("Total energy")
    println("------------")

    Ham.energies = calc_energies( Ham, psiks )
    if use_smearing
        Ham.energies.mTS = Entropy
    end
    println( Ham.energies )

end