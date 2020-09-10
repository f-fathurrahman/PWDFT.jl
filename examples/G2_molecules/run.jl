using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures", "DATA_G2_mols")

import InteractiveUtils
InteractiveUtils.versioninfo()
import Dates
println("Now = ", Dates.now())

function get_default_psp(atoms::Atoms)
    Nspecies = atoms.Nspecies
    pspfiles = Array{String}(undef,Nspecies)
    SpeciesSymbols = atoms.SpeciesSymbols
    for isp = 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        pspfiles[isp] = joinpath(DIR_PSP, PWDFT.ALL_PADE_PSP[atsymb][1])
    end
    return pspfiles
end

function do_calc( ; molname="H2O", method="SCF" )

    Random.seed!(1234)

    # Atoms
    filename = molname*".xyz"
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, filename) )

    # Initialize Hamiltonian
    pspfiles = get_default_psp(atoms)
    ecutwfc = 20.0
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

function main()
    Nargs = length(ARGS)
    if Nargs == 0
        molname = "H2O"
    else
        molname = ARGS[1]
    end
    @time do_calc(method="SCF")
    @time do_calc(method="Emin")
    
    @time do_calc(method="SCF")
    @time do_calc(method="Emin")
end

main()