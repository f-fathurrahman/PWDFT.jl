function time_calc_rhoe()

    @printf("\n")
    @printf("----------------------------------------------------------\n")
    @printf("Timing calc_rhoe (Nkpt=fcc-8x8x8) GBRV_LDA pseudopotential\n")
    @printf("----------------------------------------------------------\n")
    @printf("\n")

    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    
    pspfile_Pt = joinpath(DIR_PWDFT, "pseudopotentials", "GBRV_LDA", "pt_lda_v1.4.uspp.F.UPF")
    pspots = [PsPot_UPF(pspfile_Pt)]
    #
    ecutwfc = 15.0
    #
    options = HamiltonianOptions()
    options.dual = 6.0
    options.meshk = [8,8,8]
    options.extra_states = 4
    #
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options)

    println("Nsyms = ", Ham.sym_info.Nsyms)
    println("Nrots = ", Ham.sym_info.Nrots)

    # Shortcuts
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    CellVolume = pw.CellVolume
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Random.seed!(4321)
    psiks = rand_BlochWavefunc( Ham )

    Rhoe = zeros(Float64, Npoints, Nspin)

    @printf("Using calc_rhoe:  ")
    #@btime $Rhoe = calc_rhoe( $Ham, $psiks )
    @time Rhoe = calc_rhoe( Ham, psiks )
    println("integ Rhoe = ", sum(Rhoe)*CellVolume/Npoints)

    @printf("Using calc_rhoe!: ")
    #@btime calc_rhoe!( $Ham, $psiks, $Rhoe )
    @time calc_rhoe!( Ham, psiks, Rhoe )
    println("integ Rhoe = ", sum(Rhoe)*CellVolume/Npoints)

end

time_calc_rhoe()
