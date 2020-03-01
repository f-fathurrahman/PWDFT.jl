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
function KS_solve_SCF_NLsolve!( Ham::Hamiltonian )

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
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    evals = zeros( Float64, Nstates, Nkspin )

    ethr_last = 1e-5
    betamix = 1.0 # use Rhoe_out directly

    function density_map( Rhoe_in )
    
        update!( Ham, Rhoe_in )
    
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr_last,
                      Nstates_conv=Nstates_occ )

        Rhoe_out = calc_rhoe( Ham, psiks )

        # simple mixing here?
        Rhoe_next = betamix*Rhoe_out + (1 - betamix)*Rhoe_in
        return Rhoe_next

    end

    NiterMax = 100

    fpres = NLsolve.nlsolve(x -> density_map(x) - x, Rhoe,
            m=5, method=:anderson, xtol=1e-5, ftol=0.0, show_trace=true,
            iterations=NiterMax )
    #( fixpoint=fpres.zero, converged=NLsolve.converged(fpres) )

    #solver = scf_NLsolve_solver()
    #fpres = solver(density_map, Rhoe, 1e-5, NiterMax)
    
    update!( Ham, fpres.zero )

    println()
    println("Total energy")
    println("------------")

    Ham.energies = calc_energies( Ham, psiks )
    println( Ham.energies )

end