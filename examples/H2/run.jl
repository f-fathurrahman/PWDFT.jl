using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main( ; method="SCF" )

    Random.seed!(1234)

    # Atoms
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="pulay" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

@time main(method="SCF")
@time main(method="Emin")

@time main(method="SCF")
@time main(method="Emin")