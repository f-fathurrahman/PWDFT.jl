
function create_Ham_Pt_fcc_smearing_gbrv(; meshk=[3,3,3])
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspots = [
        PsPot_UPF(joinpath(DIR_PSP_LDA_GBRV, "pt_lda_v1.4.uspp.F.UPF"))
    ]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    #
    options = HamiltonianOptions()
    options.extra_states = 4
    options.dual = ecutrho/ecutwfc
    options.meshk = meshk
    #
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    # Set some options
    Ham.electrons.use_smearing = true
    Ham.electrons.kT = 0.003
    # Make sure that there are empty states
    @assert Ham.electrons.Nstates > Ham.electrons.Nstates_occ
    # TODO: Need to check if the states are enough?

    return Ham
end
