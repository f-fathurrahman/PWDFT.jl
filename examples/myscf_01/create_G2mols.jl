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

function init_Ham_G2mols( molname::String; xcfunc="VWN" )
    filename = molname*".xyz"
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "DATA_G2_mols", filename) )
    pspfiles = get_default_psp(atoms, xcfunc=xcfunc)
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc=xcfunc )
end
