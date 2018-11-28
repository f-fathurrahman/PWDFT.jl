function time_PWGrid()

    @printf("\n")
    @printf("----------------------\n")
    @printf("Timing PWGrid function\n")
    @printf("----------------------\n")
    @printf("\n")

    ECUTWFC = range(15.0, stop=50.0, step=5.0)    
    LatVecs = gen_lattice_sc(16.0)

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime PWGrid($ecutwfc, $LatVecs)
    end
end

time_PWGrid()