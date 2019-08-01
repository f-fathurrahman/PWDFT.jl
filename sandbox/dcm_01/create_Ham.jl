function create_Ham_Pt_fcc_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspfiles = [joinpath(DIR_PSP, "Pt-q18.gth")]
    ecutwfc = 30.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[8,8,8], extra_states=4 )
end


function create_Ham_O2_smearing()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "O2.xyz"),
                   LatVecs = gen_lattice_cubic(16.0) )
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_NH3()
    atoms = Atoms(ext_xyz_file="NH3.xyz")
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end


function create_Ham_N2()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "N2.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end


function create_Ham_CO()
    # Atoms
    atoms = Atoms(xyz_string=
        """
        2

        O    0.000000  -0.073613     0.000000
        C    0.000000   1.073613     0.000000
        """, LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth"),
                joinpath(DIR_PSP, "C-q4.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc )
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



function create_Ham_GaAs_v1()

    LatVecs = zeros(3,3)
    LatVecs[:,1] = [0.5, 0.5, 0.0]
    LatVecs[:,2] = [0.5, 0.0, 0.5]
    LatVecs[:,3] = [0.0, 0.5, 0.5]
    LatVecs = LatVecs*5.6537*ANG2BOHR

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.1  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=LatVecs)

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function create_Ham_GaAs_v2()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.1  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_Fe_bcc()
    
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, in_bohr=true, LatVecs=gen_lattice_bcc(2.856*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Fe-q8.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3],
                        Nspin=2, extra_states=4 )
end