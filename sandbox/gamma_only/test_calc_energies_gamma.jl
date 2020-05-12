using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
include("BlochWavefuncGamma.jl")
include("calc_rhoe_gamma.jl")
include("Poisson_solve_gamma.jl")
include("op_K_gamma.jl")
include("op_V_loc_gamma.jl")
include("op_V_Ps_nloc_gamma.jl")
include("op_H_gamma.jl")
include("calc_energies_gamma.jl")

include("unfold_BlochWavefuncGamma.jl")

function test_01()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )

    psis = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis)

    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    Rhoe = calc_rhoe(Ham, psis)

    Rhoe_ = calc_rhoe(Ham_, psiks)

    update!(Ham, Rhoe)
    update!(Ham_, Rhoe_)

    Ekin = calc_E_kin(Ham, psis)
    Ekin_ = calc_E_kin(Ham_, psiks)
    println("Ekin  = ", Ekin)
    println("Ekin_ = ", Ekin_)

    E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham )
    E_Ps_loc_, E_Hartree_, E_xc_ = calc_E_local( Ham_ )

    println("E_Ps_loc = ", E_Ps_loc)
    println("E_Ps_loc = ", E_Ps_loc)

    println("E_Hartree  = ", E_Hartree)
    println("E_Hartree_ = ", E_Hartree_)

    println("E_xc  = ", E_xc)
    println("E_xc_ = ", E_xc_)

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc  = calc_E_Ps_nloc( Ham, psis )
        E_Ps_nloc_ = calc_E_Ps_nloc( Ham_, psiks )
    else
        E_Ps_nloc  = 0.0
        E_Ps_nloc_ = 0.0
    end
    println("E_Ps_nloc  = ", E_Ps_nloc)
    println("E_Ps_nloc_ = ", E_Ps_nloc_)

    energies  = calc_energies(Ham, psis)
    energies_ = calc_energies(Ham_, psiks)

    println()
    println(energies)

    println()
    println(energies_)

end

test_01()