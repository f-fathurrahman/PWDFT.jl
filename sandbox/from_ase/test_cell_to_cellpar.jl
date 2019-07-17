using PWDFT
using LinearAlgebra

include("cell_to_cellpar.jl")

function test_main()
    #LatVecs = gen_lattice_fcc(10.0)
    #LatVecs = gen_lattice_sc(10.0)
    #LatVecs = gen_lattice_bcc(10.0)
    #LatVecs = gen_lattice_bcc_v2(10.0)
    #LatVecs = gen_lattice_rhombohedral(10.0, 30.0)
    LatVecs = gen_lattice_triclinic(10.0, 12.0, 11.0, 80.0, 40.0, 70.0)


    cell_to_cellpar(LatVecs)
end

test_main()