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
include("unfold_BlochWavefuncGamma.jl")

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

    psis = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis)

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    Rhoe = calc_rhoe(Ham, psis, renormalize=false)
    integRhoe = sum(Rhoe)*dVol
    println("integRhoe (with renormalize=false) = ", integRhoe)

    Rhoe = calc_rhoe(Ham, psis, renormalize=true)
    integRhoe = sum(Rhoe)*dVol
    println("integRhoe (with renormalize=true) = ", integRhoe)

    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
    PWDFT.ortho_check(psiks[1])

    Rhoe_ = calc_rhoe(Ham_, psiks)

    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Rhoe[ip], Rhoe_[ip])
    end

end

test_01()