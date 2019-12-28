using Random

using CuArrays
using PWDFT
using PWDFT_cuda

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main_CPU()

    Random.seed!(1234)

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    Nspin = 1
    
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    println( Ham_cpu )

    KS_solve_Emin_PCG!( Ham_cpu, skip_initial_diag=true, startingrhoe=:random )
    
    @time KS_solve_Emin_PCG!( Ham_cpu, skip_initial_diag=true, startingrhoe=:random )

end

function main_GPU()

    Random.seed!(1234)

    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    Nspin = 1

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    KS_solve_Emin_PCG!( Ham, skip_initial_diag=true, startingrhoe=:random )
    
    @time KS_solve_Emin_PCG!( Ham, skip_initial_diag=true, startingrhoe=:random )

end


main_CPU()
main_GPU()
