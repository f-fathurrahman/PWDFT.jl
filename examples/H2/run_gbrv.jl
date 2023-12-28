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
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "GBRV_LDA", "h_lda_v1.4.uspp.F.UPF")]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, dual=dual )
    println(Ham)

    # XXX: This should take into account the overlap operator
    # XXX: Currently this is handled in _prepare_scf (called electrons_scf)
    psiks = rand_BlochWavefunc(Ham)

    # XXX: Only this solver is currently implemented
    electrons_scf!(Ham, psiks)
end

@time main()
@time main()
