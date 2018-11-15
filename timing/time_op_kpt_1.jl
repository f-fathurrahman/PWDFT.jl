function time_op_kpt_1()

    @printf("\n")
    @printf("-------------------------------------\n")
    @printf("Timing Hamiltonian operators (Nkpt=1)\n")
    @printf("-------------------------------------\n")
    @printf("\n")

    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[8,8,8], extra_states=4 )
        
    # Shortcuts
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    CellVolume = pw.CellVolume
    Ns = pw.Ns

    Random.seed!(4321)
    psiks = rand_BlochWavefunc( Ham )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    @printf("Integ rhoe = %18.10f\n\n", sum(Rhoe)*CellVolume/prod(Ns))

    @printf("Calc rhoe                     : ")
    @btime calc_rhoe( $Ham, $psiks )

    @printf("Updating local potential      : ")
    @btime update!( $Ham, $Rhoe )

    @printf("Kinetic operator              : ")
    @btime op_K( $Ham, $psiks )

    @printf("V_Ps_loc                      : ")
    @btime op_V_Ps_loc( $Ham, $psiks )

    @printf("V_loc (Ps loc + Hartree + XC) : ")
    @btime op_V_loc( $Ham, $psiks )

    @printf("V_Ps_nloc                     : ")
    @btime op_V_Ps_nloc( $Ham, $psiks )

    @printf("Hamiltonian                   : ")
    @btime op_H( $Ham, $psiks )

end