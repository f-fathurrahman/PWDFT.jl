using Printf
using PWDFT
using LinearAlgebra
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("my_scf.jl")
include("create_G2mols.jl")

function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    #pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    pspfiles = [ "Si.upf" ]
    ecutwfc = 15.0
    #return Hamiltonian( atoms, pspfiles, ecutwfc,
    #    meshk=[3,3,3], Ns_=(32,32,32) )
    return Hamiltonian( atoms, pspfiles, ecutwfc,
        meshk=[3,3,3] )
end

function main()
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    #Ham = init_Ham_G2mols("H2O")
    println(Ham)
    psiks = rand_BlochWavefunc(Ham)
    my_scf!(Ham, psiks, betamix=0.5)
    #KS_solve_SCF!(Ham, psiks, mix_method="broyden")
    #KS_solve_Emin_PCG!(Ham, psiks)
end

main()

