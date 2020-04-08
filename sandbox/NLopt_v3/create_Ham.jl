function create_Ham_ZnO()
    DIR_PWINPUT = joinpath(DIR_STRUCTURES, "DATA_Semiconductors")
    atoms, meshk = read_pwscf_input(joinpath(DIR_PWINPUT, "ZnO.in"))

#    atoms = Atoms( xyz_string_frac=
#        """
#        4
#
#        Zn      0.3333333   0.6666667   0.0000000
#        Zn      0.6666667   0.3333333   0.5000000
#        O       0.3333333   0.6666667   0.3450000
#        O       0.6666667   0.3333333   0.8450000
#        """, in_bohr=true,
#        LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR ) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Zn-q2.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    #return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[1,1,1] )
end

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
            #return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], use_symmetry=false, time_reversal=false )
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
        end
    end
end

function create_Ham_GaAs()

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(10.6839444516) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    #return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], use_symmetry=false, time_reversal=false )
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_GaAs_cubic()
    DIR_PWINPUT = joinpath(DIR_STRUCTURES, "DATA_Semiconductors")
    atoms, meshk = read_pwscf_input("GaAs.in")
    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "As-q5.gth"),
                joinpath(DIR_PSP, "Ga-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
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

function create_Ham_NH3()
    # Atoms
    atoms = init_atoms_xyz(joinpath(DIR_STRUCTURES,"NH3.xyz"))
    atoms.LatVecs = gen_lattice_cubic(16.0)

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    return Ham
end
