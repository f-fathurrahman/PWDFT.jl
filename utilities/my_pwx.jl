using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function my_pwx(; filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

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

    electrons_scf!(Ham, psiks, NiterMax=100, use_smearing=use_smearing, kT=kT, betamix=0.1)

#=
    # Not yet working for smearing
    KS_solve_Emin_PCG!(Ham, psiks, NiterMax=100)
    energies = Ham.energies
    println()
    println(">>>> Final result:")
    println()
    
    println("-------------------------------------")
    println("Energy components in Ry")
    println("-------------------------------------")
    
    @printf("Kinetic    energy: %18.10f Ry\n", 2*energies.Kinetic )
    @printf("Ps_loc     energy: %18.10f Ry\n", 2*energies.Ps_loc )
    @printf("Ps_nloc    energy: %18.10f Ry\n", 2*energies.Ps_nloc )
    @printf("Hartree    energy: %18.10f Ry\n", 2*energies.Hartree )
    @printf("XC         energy: %18.10f Ry\n", 2*energies.XC )
    @printf("-TS              : %18.10f Ry\n", 2*energies.mTS)
    @printf("-------------------------------------\n")
    
    E_elec = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
             energies.Hartree + energies.XC + energies.mTS
    
    @printf("Electronic energy: %18.10f Ry\n", 2*E_elec)
    @printf("NN         energy: %18.10f Ry\n", 2*energies.NN )
    @printf("-------------------------------------\n")

    E_total = E_elec + energies.NN
    @printf("! Total = %18.10f Ry\n", 2*E_total)
=#

end

my_pwx()
