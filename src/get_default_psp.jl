# Return Vector of PsPot_GTH for given `atoms`
function get_default_PsPot_GTH( atoms::Atoms; xcfunc="VWN" )
    DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
    if xcfunc == "VWN"
        DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
        ALL_PSP = PWDFT.ALL_PADE_PSP
    elseif xcfunc == "PBE"
        DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth")
        ALL_PSP = PWDFT.ALL_PBE_PSP
    else
        error("xcfunc=$xcfunc in get_default_psp is not known")
    end
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    pspots = Vector{PsPot_GTH}(undef,Nspecies)
    for isp in 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        pspots[isp] = PsPot_GTH(joinpath(DIR_PSP, ALL_PSP[atsymb][1]))
    end
    return pspots
end