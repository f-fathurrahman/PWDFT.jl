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

        # The result is far better for v.1.4.0-rc2 as compared to v1.3.1
        @printf("find_symm_bravais_latt: ")
        @btime Nrots, s, sname = PWDFT.find_symm_bravais_latt( $atoms.LatVecs )
        println()
    end
end

time_SymmetryInfo()

