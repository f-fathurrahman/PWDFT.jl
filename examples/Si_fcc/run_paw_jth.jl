using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

import InteractiveUtils
InteractiveUtils.versioninfo()
import Dates
println("Now = ", Dates.now())

function main( ; method="SCF" )

    Random.seed!(1234)

    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "PAW_JTH_LDA", "Si.upf")]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, dual=dual, meshk=[3,3,3] )
    println(Ham)

    psiks = rand_BlochWavefunc(Ham)

    # XXX: Only this solver is currently implemented
    electrons_scf!(Ham, psiks=psiks)
end

@time main()
@time main()
