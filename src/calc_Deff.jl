function _calc_Deff!(
    ispin::Int64, atoms, pspotNL, ebnd::Float64, Deff
)
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    Deeq = pspotNL.Deeq
    qq_at = pspotNL.qq_at
    #
    ok_uspp_or_paw = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)
    #
    @views Deff[:,:,:] .= Deeq[:,:,:,ispin]
    #println("ok_uspp_or_paw = ", ok_uspp_or_paw)
    #
    if ok_uspp_or_paw
        for isp in 1:Nspecies, ia in 1:Natoms
            if isp != atm2species[ia]
                continue
            end
            Deff[:,:,ia] .-= ebnd * qq_at[:,:,ia]
        end
    end
    #println("sum qq_at = ", sum(qq_at))
    #println("sum Deeq = ", sum(Deeq))
    #println("sum Deff = ", sum(Deff))
    return
end
