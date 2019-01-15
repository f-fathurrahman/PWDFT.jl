function time_Hamiltonian_H()
    
    @printf("\n")
    @printf("--------------------------------------------\n")
    @printf("Timing Hamiltonian initialization for atom H\n")
    @printf("--------------------------------------------\n")
    @printf("\n")

    atoms = Atoms( xyz_string="""
        1
    
        H   0.0   0.0   0.0
        """,
        LatVecs=gen_lattice_sc(16.0) )
    
    pspfiles = [joinpath(DIR_PSP,"H-q1.gth")]
    ECUTWFC = range(15.0, stop=50.0, step=5.0)    

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime Hamiltonian( $atoms, $pspfiles, $ecutwfc )
    end

end

time_Hamiltonian_H()