function time_gen_lattice()

    @printf("\n")
    @printf("------------------------------\n")
    @printf("Timing gen_lattice_* functions\n")
    @printf("------------------------------\n")
    @printf("\n")

    @printf("gen_lattice_sc  : ")
    @btime gen_lattice_sc(10.0)

    @printf("gen_lattice_fcc : ")
    @btime gen_lattice_fcc(10.0)
        
    @printf("gen_lattice_bcc : ")
    @btime gen_lattice_bcc(10.0)
end