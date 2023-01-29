using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

include("PWSCFInput.jl")
include("init_Ham_from_pwinput.jl")

function test_main()
    Ham, pwinput = init_Ham_from_pwinput()

    write_xsf("ATOMS_from_pwinput.xsf", Ham.atoms)
    println(Ham)

    # XXX: This should take into account the overlap operator
    # XXX: Currently this is handled in _prepare_scf (called electrons_scf)
    psiks = rand_BlochWavefunc(Ham)

    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
    end

    electrons_scf!(Ham, psiks, NiterMax=100, use_smearing=use_smearing, kT=kT)
end

test_main()
