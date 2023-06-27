#=
Collection of functions that create some objects that might be used
in several test.
=#


function create_atoms_Si_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    return atoms
end


function create_Ham_Si_fcc_oncv()
    atoms = create_atoms_Si_fcc()

    # Initialize Hamiltonian
    pspots = [
        PsPot_UPF(joinpath(DIR_PWDFT, "pseudopotentials",
            "ONCV_v0.4.1_LDA", "Si.upf"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    options = HamiltonianOptions()
    options.meshk = [3,3,3]
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end


function create_Ham_Si_fcc_gbrv()
    atoms = create_atoms_Si_fcc()

    # Initialize Hamiltonian
    pspots = [
        PsPot_UPF(joinpath(DIR_PWDFT, "pseudopotentials",
            "GBRV_LDA", "si_lda_v1.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    options.meshk = [3,3,3]
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end


function create_Ham_Si_fcc_paw_pslib()
    atoms = create_atoms_Si_fcc()

    # Initialize Hamiltonian
    pspots = [
        PsPot_UPF(joinpath(DIR_PWDFT, "pseudopotentials",
            "PSLIB_US_PAW_LDA", "Si.pz-n-kjpaw_psl.1.0.0.UPF"))
    ]
    # FIXME: PSLIB_US_PAW_LDA is not included in the repo
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    options.meshk = [3,3,3]
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end

