import InteractiveUtils
using Printf
using BenchmarkTools
using PWDFT

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


function time_PWGrid()

    @printf("\n")
    @printf("-----------------------\n")
    @printf("Timing PWGrid functions\n")
    @printf("-----------------------\n")
    @printf("\n")

    ECUTWFC = range(15.0, stop=50.0, step=5.0)    
    LatVecs = gen_lattice_sc(16.0)

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime PWGrid($ecutwfc, $LatVecs)
    end
end

println()
InteractiveUtils.versioninfo()

time_gen_lattice()
time_Atoms()
time_PWGrid()

