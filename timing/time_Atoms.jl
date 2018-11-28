function time_Atoms()

    @printf("\n")
    @printf("----------------------\n")
    @printf("Timing Atoms functions\n")
    @printf("----------------------\n")
    @printf("\n")

    @printf("Using default constructor (dummy Atoms)          : ")
    @btime Atoms()

    @printf("Using xyz_string with two atoms, default LatVecs : ")
    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """)

    @printf("Using xyz_string with two atoms, given LatVecs   : ")
    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """, LatVecs=gen_lattice_sc(10.0))
end

time_Atoms()