using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

import NLsolve


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


function create_Ham_Si_fcc()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    return Ham
end


function main_scf()
    Random.seed!(1234)
    Ham = create_Ham_Si_fcc()
    KS_solve_SCF!( Ham, mix_method="simple", betamix=0.5 )
end


function main()

    Random.seed!(1234)

    Ham = create_Ham_Si_fcc()

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
    betamix = 1.0

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
    solver = scf_NLsolve_solver()
    fpres = solver(density_map, Rhoe, 1e-7, NiterMax)

    println(fpres.converged)
    println(size(fpres.fixpoint))
    println(Npoints)
    
    update!( Ham, fpres.fixpoint )
    Ham.energies = calc_energies( Ham, psiks )
    println( Ham.energies )

    println("Pass here")
end

main()

#main_scf()