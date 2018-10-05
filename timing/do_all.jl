using BenchmarkTools
using PWDFT

function time_01()

    @btime Atoms()


    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """)


    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """, LatVecs=gen_lattice_sc(10.0))

    @btime gen_lattice_sc(10.0)
    @btime gen_lattice_fcc(10.0)
    @btime gen_lattice_bcc(10.0)
end


function time_02()

    ECUTWFC = ( 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0 )    
    LatVecs = gen_lattice_sc(16.0)

    for ecutwfc in ECUTWFC
        @btime PWGrid($ecutwfc, $LatVecs)
    end
end


#time_01()
time_02()

