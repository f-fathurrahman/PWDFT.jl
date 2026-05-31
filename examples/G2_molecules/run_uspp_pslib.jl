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
        dir_psp = joinpath(DIR_PWDFT, "pseudopotentials", "PSLIB_US_PAW_PBE")
    elseif xcfunc == "VWN"
        dir_psp = joinpath(DIR_PWDFT, "pseudopotentials", "PSLIB_US_PAW_LDA")
    else
        error("Unsupported xcfunc = $xcfunc")
    end

    #XXX These can be cached ?
    listfiles = readdir(dir_psp)
    atsymb_list = String[]
    for l in listfiles
        atsymb = split(l, ".")[1]
        if !(atsymb in atsymb_list)
            append!(atsymb_list, [atsymb])
        end
    end
    dict_kjpaw = Dict{String,Vector{String}}()
    for atsymb in atsymb_list
        list_psp_atsymb = String[]
        for l in listfiles
            if occursin("rrkjus", l) && atsymb==split(l, ".")[1]
                append!(list_psp_atsymb, [l])
            end
        end
        merge!(dict_kjpaw, Dict(atsymb => list_psp_atsymb))
    end

    Nspecies = atoms.Nspecies
    pspfiles = Array{String}(undef,Nspecies)
    SpeciesSymbols = atoms.SpeciesSymbols
    for isp in 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        # choose the one with the shortest filename
        idx_chosen = argmin(length.(dict_kjpaw[atsymb]))
        pspfiles[isp] = joinpath(dir_psp, dict_kjpaw[atsymb][idx_chosen])
    end
    return pspfiles
end



function do_calc( molname::String; xcfunc="VWN" )

    Random.seed!(1234)

    # Atoms
    filename = molname*".xyz"
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, filename) )

    # Initialize Hamiltonian
    pspfiles = get_default_psp(atoms, xcfunc=xcfunc)
    ecutwfc = 20.0
    dual = 5.0
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc, dual = dual, xcfunc = xcfunc
    )
    println(Ham)
    #
    psiks = rand_BlochWavefunc(Ham)
    electrons_scf!(Ham, psiks=psiks)
    return
end

function main()
    Nargs = length(ARGS)
    if Nargs == 0
        molname = "H2O"
    else
        molname = ARGS[1]
    end
    do_calc(molname, xcfunc="VWN")
end

main()