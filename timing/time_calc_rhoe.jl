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
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin

    Random.seed!(4321)
    psiks = rand_BlochWavefunc( Ham )

    Rhoe = zeros(Npoints, Nspin)

    @printf("Using calc_rhoe:  ")
    @btime $Rhoe = calc_rhoe( $Ham, $psiks )

    @printf("Using calc_rhoe!: ")
    @btime calc_rhoe!( $Ham, $psiks, $Rhoe )
end


function time_calc_rhoe2()

    @printf("\n")
    @printf("------------------------------\n")
    @printf("Timing calc_rhoe (CO molecule)\n")
    @printf("------------------------------\n")
    @printf("\n")

    atoms = Atoms(xyz_string_frac=
        """
        2

        C  0.0  0.0  0.0
        O  1.5  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    
    pspfiles = [joinpath(DIR_PSP, "C-q4.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
        
    # Shortcuts
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    CellVolume = pw.CellVolume
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin

    Random.seed!(4321)
    psiks = rand_BlochWavefunc( Ham )

    Rhoe = zeros(Npoints,Nspin)

    @printf("Using calc_rhoe:  ")
    @btime $Rhoe = calc_rhoe( $Ham, $psiks )

    @printf("Using calc_rhoe!: ")
    @btime calc_rhoe!( $Ham, $psiks, $Rhoe )

end

time_calc_rhoe()
time_calc_rhoe2()

