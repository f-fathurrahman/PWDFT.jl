function create_Ham_Si_fcc( ; xcfunc="VWN", Nspin=1 )

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]

    ecutwfc = 15.0
    if xcfunc == "PBE"
        if Nspin == 2
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], xcfunc="PBE" )
        else
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], xcfunc="PBE", Nspin=2, extra_states=4 )
        end
    else
        if Nspin == 2
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2, extra_states=4 )
        else
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
        end
    end
end

function create_Ham_H2()
    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function create_Ham_H_atom()
    atoms = Atoms(xyz_string=
        """
        1

        H      0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end