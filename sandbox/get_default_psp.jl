function get_default_psp(atoms::Atoms)
    DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
    DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

    Nspecies = atoms.Nspecies
    pspfiles = Array{String}(undef,Nspecies)

    SpeciesSymbols = atoms.SpeciesSymbols
    for isp = 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        pspfiles[isp] = joinpath(DIR_PSP, PWDFT.ALL_PADE_PSP[atsymb][1])
    end
    return pspfiles
end