function time_Hamiltonian()
    
    @printf("\n")
    @printf("---------------------------\n")
    @printf("Timing Hamiltonian function\n")
    @printf("---------------------------\n")
    @printf("\n")

    atoms = Atoms( xyz_string="""
        1
    
        H   0.0   0.0   0.0
        """,
        LatVecs=gen_lattice_sc(16.0) )
    
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]    
    ECUTWFC = range(15.0, stop=50.0, step=5.0)    

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime Hamiltonian( $atoms, $pspfiles, $ecutwfc )
    end

end