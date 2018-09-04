function calc_PspCore_ene( atoms::Atoms, PsPots::Array{PsPot_GTH,1}, CellVolume::Float64 )

    Nspecies = atoms.Nspecies
    epsatm = zeros(Nspecies)

    for isp = 1:Nspecies
        zval = atoms.Zvals[isp]
        psp = PsPots[isp]
        rloc = psp.rlocal
        zval = psp.zval
        c1 = psp.c[1]
        c2 = psp.c[2]
        c3 = psp.c[3]
        c4 = psp.c[4]
        epsatm[isp] = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15.0*c3 + 105.0*c4)
    end

    chg = 0.0
    Natoms = atoms.Natoms
    pspcore_ene = 0.0
    for ia = 1:Natoms
        isp = atoms.atm2species[ia]
        chg = chg + atoms.Zvals[isp]
        pspcore_ene = pspcore_ene + epsatm[isp]
    end
    
    pspcore_ene = chg*pspcore_ene/CellVolume
    return pspcore_ene

end