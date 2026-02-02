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
        PsPot_UPF(joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Si.upf"))
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
        PsPot_UPF(joinpath(DIR_PSP, "GBRV_LDA", "si_lda_v1.uspp.F.UPF"))
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


function create_atoms_N2H4()
    atoms = Atoms(xyz_string="""
6

N       5.94821400       6.81171100       5.22639100
N       5.94821400       5.37379300       5.22639100
H       6.15929600       7.18550400       6.15196500
H       5.00000000       7.09777800       5.00000000
H       5.73713200       5.00000000       6.15196500
H       6.89642800       5.08772600       5.00000000
""", LatVecs=gen_lattice_sc(16.0))
    return atoms
end

function create_Ham_N2H4_oncv()
    atoms = create_atoms_N2H4()
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "N.upf")),
        PsPot_UPF(joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "H.upf"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    options = HamiltonianOptions()
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end



function create_Ham_N2H4_gbrv()
    atoms = create_atoms_N2H4()
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP, "GBRV_LDA", "n_lda_v1.2.uspp.F.UPF")),
        PsPot_UPF(joinpath(DIR_PSP, "GBRV_LDA", "h_lda_v1.4.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    options = HamiltonianOptions()
    options.dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    return Ham
end


function create_Ham_Fe_bcc_oncv()
    atoms = Atoms(2, 1,
        [0.0 2.7079775377810678; 0.0 2.7079775377810678; 0.0 2.7079775377810678],
        [1, 1], ["Fe", "Fe"], ["Fe"],
        [5.4159550755621355 0.0 0.0; 0.0 5.4159550755621355 0.0; 0.0 0.0 5.4159550755621355],
        [16.0], [0.0]
    );
    pspfiles = [joinpath(DIR_PSP, "ONCV_v0.4.1_LDA", "Fe.upf")];
    ecutwfc = 20.0;
    options = HamiltonianOptions(4.0,
        2, 2, [3, 3, 3], [0, 0, 0], true, (0, 0, 0), nothing, "", "VWN",
        false, -1, 20, true, true, 0.001, [0.4], nothing, nothing, false, false
    );
    #(;(v=>getfield(atoms, v) for v in fieldnames(typeof(atoms)))...)
    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies);
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp]);
    end
    Ham = Hamiltonian(atoms, pspots, ecutwfc, options);
    return Ham
end