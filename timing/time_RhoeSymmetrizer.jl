function time_RhoeSymmetrizer()
    @printf("\n")
    @printf("-------------------------------\n")
    @printf("Timing RhoeSymmetrizer function\n")
    @printf("-------------------------------\n")
    @printf("\n")
    
    DIR_PWINPUT = joinpath(DIR_STRUCTURES, "DATA_DeltaCodes", "PWINPUT")
    ecutwfc = 15.0

    for fil in ["Si.in", "Pt.in", "N.in", "H.in"]

        atoms, _ = read_pwscf_input(joinpath(DIR_PWINPUT, fil))
        pw = PWGrid(ecutwfc, atoms.LatVecs)
        sym_info = SymmetryInfo(atoms)

        @printf("File = %5s, Natoms = %d, Nsyms = %d : ", fil, atoms.Natoms, sym_info.Nsyms )
        @btime rhoe_sym = RhoeSymmetrizer($atoms, $pw, $sym_info)

    end
end

time_RhoeSymmetrizer()
