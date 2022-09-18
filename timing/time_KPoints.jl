function time_KPoints()
    @printf("\n")
    @printf("-----------------------\n")
    @printf("Timing KPoints function\n")
    @printf("-----------------------\n")
    @printf("\n")
    
    DIR_PWINPUT = joinpath(DIR_STRUCTURES, "DATA_DeltaCodes", "PWINPUT")
    atoms, meshk = read_pwscf_input(joinpath(DIR_PWINPUT, "Si.in"))

    sym_info = SymmetryInfo(atoms)
    shiftk = [0,0,0]

    println("Internal (pure Julia): ")
    for k in range(5,stop=10)
        meshk = [k, k, k]
        @printf("meshk = [%2d,%2d,%2d] : ", k, k, k)
        @btime kpoints = KPoints( $atoms, $meshk, $shiftk, $sym_info.s )
    end

end

time_KPoints()