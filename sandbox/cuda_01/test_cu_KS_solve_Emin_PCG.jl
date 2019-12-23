using Test

include("PWDFT_cuda.jl")

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main()

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 50.0

    Nspin = 1

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    
    Random.seed!(1234)
    KS_solve_Emin_PCG!( Ham, skip_initial_diag=true, startingrhoe=:random )

    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    println( Ham_cpu )
    
    Random.seed!(1234)
    KS_solve_Emin_PCG!( Ham_cpu, skip_initial_diag=true, startingrhoe=:random )

    Random.seed!(1234)
    @time KS_solve_Emin_PCG!( Ham, skip_initial_diag=true, startingrhoe=:random )
    
    Random.seed!(1234)
    @time KS_solve_Emin_PCG!( Ham_cpu, skip_initial_diag=true, startingrhoe=:random )

end

main()
