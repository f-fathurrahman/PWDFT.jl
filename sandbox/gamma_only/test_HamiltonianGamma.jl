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

function test_01()
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

    println("sum V Ps loc gamma = ", sum(Ham.potentials.Ps_loc))
    println("sum V Ps loc       = ", sum(Ham_.potentials.Ps_loc))

    println("betaNL gamma 1 = ", Ham.pspotNL.betaNL[1,1])
    println("betaNL       1 = ", sum(Ham_.pspotNL.betaNL[1][1,1]))

    println(size(Ham.pspotNL.betaNL))
    println(size(Ham_.pspotNL.betaNL[1]))

    s1 = sum(Ham.pspotNL.betaNL[2:end])
    s2 = sum(Ham_.pspotNL.betaNL[1][2:end])
    println("s1 = ", s1)
    println("s2 = ", s2)
    println("s2/s1 = ", s2/s1)

    betaNL = Ham.pspotNL.betaNL
    println( 2*dot(betaNL,betaNL) - conj(betaNL[1])*betaNL[1] )

    betaNL1 = Ham_.pspotNL.betaNL[1]
    println(dot(betaNL1,betaNL1))

    psis = randn_BlochWavefuncGamma(Ham)
    Rhoe = calc_rhoe(Ham, psis)

    update!(Ham, Rhoe)

end

test_01()