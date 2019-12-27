using Test

include("PWDFT_cuda.jl")

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main_CPU()

    Random.seed!(1234)

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "PtO.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    pspfiles = [joinpath(DIR_PSP, "Pt-q18.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 50.0

    Nspin = 1

    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    println( Ham_cpu )
    
    Random.seed!(1234)
    KS_solve_Emin_PCG!( Ham_cpu, skip_initial_diag=true, startingrhoe=:random )

end

main_CPU()
