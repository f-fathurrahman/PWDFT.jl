using Printf
using PWDFT

function test_main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(10.2631)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3] )

    psp = Ham.pspots[1]
    rloc = psp.rlocal
    zval = psp.zval
    c1 = psp.c[1]
    c2 = psp.c[2]
    c3 = psp.c[3]
    c4 = psp.c[4]
    epsatm = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15.0*c3 + 105.0*c4)
    chg = 0.0
    Natoms = atoms.Natoms
    for ia = 1:Natoms
        isp = atoms.atm2species[ia]
        chg = chg + atoms.Zvals[isp]
    end
    pspcore_ene = epsatm*Natoms*chg/Ham.pw.CellVolume

    println("pspcore_ene = ", pspcore_ene)


end

@time test_main()
