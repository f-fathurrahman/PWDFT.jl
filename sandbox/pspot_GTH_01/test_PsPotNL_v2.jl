using LinearAlgebra
using Random
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("PsPotNL_v2.jl")

function test_main()

    atoms = Atoms(xyz_string="""
        1

        Cs  10.0  10.0  10.0
        """, in_bohr=true,
        LatVecs=gen_lattice_sc(20.0))
    println(atoms)

    ecutwfc = 80.0
    pw = PWGrid(ecutwfc, atoms.LatVecs)
    println(pw)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    Nspecies = atoms.Nspecies
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    # !!! HARDCODED !!!
    pspots[1] = PsPot_GTH(joinpath(DIR_PSP,"Cs-q9.gth"))

    for psp in pspots
        println(psp)
    end

    pspotNL = PsPotNL_v2( atoms, pw, pspots, check_norm=true )

end

test_main()
