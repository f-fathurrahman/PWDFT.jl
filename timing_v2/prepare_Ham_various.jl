function create_Ham_atom_Pt_gth()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspots = [ PsPot_GTH(joinpath(DIR_PSP_LDA_GTH, "Pt-q18.gth")) ]
    ecutwfc = 20.0
    options = HamiltonianOptions()
    return Hamiltonian(atoms, pspots, ecutwfc, options)
end

function create_Ham_atom_Pt_oncv()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspots = [ PsPot_UPF(joinpath(DIR_PSP_LDA_ONCV, "Pt.upf")) ]
    ecutwfc = 20.0
    options = HamiltonianOptions()
    return Hamiltonian(atoms, pspots, ecutwfc, options)
end

function create_Ham_atom_Pt_gbrv()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspots = [ PsPot_UPF(joinpath(DIR_PSP_LDA_GBRV, "pt_lda_v1.4.uspp.F.UPF")) ]
    ecutwfc = 20.0
    ecutrho = 100.0
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    return Hamiltonian(atoms, pspots, ecutwfc, options)
end


function create_Ham_atom_Pt_paw_jth()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspots = [ PsPot_UPF(joinpath(DIR_PSP_LDA_PAW_JTH, "Pt.upf")) ]
    ecutwfc = 20.0
    ecutrho = 100.0
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    return Hamiltonian(atoms, pspots, ecutwfc, options)
end


function create_Ham_O2_gth_spinpol()
    atoms = Atoms(ext_xyz_file=joinpath(DIR_STRUCTURES, "O2.xyz"))
    pspfiles = [joinpath(DIR_PSP_LDA_GTH, "O-q6.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4, Nspin=2 )
end

function create_Ham_Pt_fcc_smearing(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspfiles = [joinpath(DIR_PSP_LDA_GTH, "Pt-q10.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=meshk, extra_states=5 )
end


function create_Ham_Pt_fcc_smearing_gbrv(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_GBRV_LDA, "pt_lda_v1.4.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    #
    return Hamiltonian( atoms, pspots, ecutwfc, options )
end


function create_Ham_Fe_bcc_smearing_gbrv(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, LatVecs=gen_lattice_bcc(2.866*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_LDA_GBRV, "fe_lda_v1.5.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    options.Nspin = 2
    #
    return Hamiltonian( atoms, pspots, ecutwfc, options )
end


function create_Ham_Fe_bcc_smearing_paw_jth(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, LatVecs=gen_lattice_bcc(2.866*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_LDA_PAW_JTH, "Fe.upf"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    options.Nspin = 2
    #
    return Hamiltonian( atoms, pspots, ecutwfc, options )
end
