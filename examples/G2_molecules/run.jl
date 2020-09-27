using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures", "DATA_G2_mols")

import InteractiveUtils
InteractiveUtils.versioninfo()
import Dates
println("Now = ", Dates.now())

function get_default_psp(atoms::Atoms; xcfunc="VWN")

    if xcfunc == "PBE"
        DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth")
        PSP_SET = PWDFT.ALL_PBE_PSP
    elseif xcfunc == "VWN"
        DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
        PSP_SET = PWDFT.ALL_PADE_PSP
    else
        error(@sprintf("Unknown xcfunc = %s", xcfunc))
    end

    Nspecies = atoms.Nspecies
    pspfiles = Array{String}(undef,Nspecies)
    SpeciesSymbols = atoms.SpeciesSymbols
    for isp = 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        pspfiles[isp] = joinpath(DIR_PSP, PSP_SET[atsymb][1])
    end
    return pspfiles
end

function do_calc( molname::String; method="SCF", xcfunc="VWN" )

    Random.seed!(1234)

    # Atoms
    filename = molname*".xyz"
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, filename) )

    # Initialize Hamiltonian
    pspfiles = get_default_psp(atoms, xcfunc=xcfunc)
    ecutwfc = 20.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
        use_xc_internal=true, xcfunc=xcfunc )
    println(Ham)
    psiks = rand_BlochWavefunc(Ham)

    if method == "SCF"
        KS_solve_SCF!(Ham, psiks, mix_method="anderson",
            startingrhoe=:random)

    elseif method == "SCF_potmix"
        KS_solve_SCF_potmix!(Ham, psiks, startingrhoe=:random)

    elseif method == "Emin"
        KS_solve_Emin_PCG!(Ham, psiks,
            startingrhoe=:random, skip_initial_diag=true)

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

function main()
    Nargs = length(ARGS)
    if Nargs == 0
        molname = "H2O"
    else
        molname = ARGS[1]
    end
    for i in 1:2
        @time do_calc(molname, method="SCF", xcfunc="PBE")
        @time do_calc(molname, method="SCF_potmix", xcfunc="PBE")
        @time do_calc(molname, method="Emin", xcfunc="PBE")
    end
end

main()