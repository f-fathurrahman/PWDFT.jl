function time_SymmetryInfo()
    @printf("\n")
    @printf("----------------------------\n")
    @printf("Timing SymmetryInfo function\n")
    @printf("----------------------------\n")
    @printf("\n")
    
    DIR_PWINPUT = joinpath(DIR_STRUCTURES, "DATA_DeltaCodes", "PWINPUT")
    
    for fil in ["Si.in", "Pt.in", "N.in", "H.in"]
        atoms, _ = read_pwscf_input(joinpath(DIR_PWINPUT, fil))
        @printf("File = %5s, Natoms = %5d : ", fil, atoms.Natoms)
        @btime sym_info = SymmetryInfo($atoms)
    end
end

time_SymmetryInfo()