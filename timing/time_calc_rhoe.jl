function time_calc_rhoe()

    @printf("\n")
    @printf("---------------------------------\n")
    @printf("Timing calc_rhoe (Nkpt=fcc-8x8x8)\n")
    @printf("---------------------------------\n")
    @printf("\n")

    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    
    pspfiles = [joinpath(DIR_PSP, "Pt-q18.gth")]
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

    Rhoe = zeros(prod(Ham.pw.Ns),2)

    @printf("Using calc_rhoe:  ")
    @btime $Rhoe = calc_rhoe( $Ham, $psiks )

    @printf("Using calc_rhoe!: ")
    @btime calc_rhoe!( $Ham, $psiks, $Rhoe )
end

time_calc_rhoe()