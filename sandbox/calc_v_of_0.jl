function calc_v_of_0( Ham::Hamiltonian )

    atoms = Ham.atoms
    pspots = Ham.pspots

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms

    v_of_0 = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        rloc = psp.rlocal
        zval = psp.zval
        c1 = psp.c[1]
        c2 = psp.c[2]
        c3 = psp.c[3]
        c4 = psp.c[4]
        v_of_0 = v_of_0 + 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15.0*c3 + 105.0*c4)  # to match QE
    end
    v_of_0 = v_of_0 / Ham.pw.CellVolume * prod(Ham.pw.Ns) / Ham.pw.CellVolume
    return v_of_0
end