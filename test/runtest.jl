using PWDFT
using Base.Test

@test gen_lattice_cubic(16.0) == 16.0*diagm(ones(3))

function test_gen_lattice()
    LatVecs1 = gen_lattice_cubic(10.0)
    println("simple cubic")
    print_matrix(LatVecs1)

    LatVecs2 = gen_lattice_fcc(10.0)
    println("fcc")
    print_matrix(LatVecs2)

    LatVecs3 = gen_lattice_bcc(10.0)
    println("bcc")
    print_matrix(LatVecs3)
    
    LatVecs4 = gen_lattice_bcc_v2(10.0)
    println("bcc v2")
    print_matrix(LatVecs4)
    
    @assert det(LatVecs3) == det(LatVecs4)

    LatVecs5 = gen_lattice_hexagonal( 10.0, 15.0 )
    println("hexagonal")
    print_matrix(LatVecs5)

    LatVecs6 = gen_lattice_trigonal( 10.0, 90.0 )
    println("trigonal")
    print_matrix(LatVecs6)

    LatVecs7 = gen_lattice_trigonal_v2( 10.0, 90.0 )
    println("trigonal v2")
    print_matrix(LatVecs7)

    LatVecs8 = gen_lattice_tetragonal_P( 10.0, 5.0 )
    println("tetragonal")
    print_matrix(LatVecs8)

    LatVecs9 = gen_lattice_tetragonal_I( 10.0, 5.0 )
    println("tetragonal v2")
    print_matrix(LatVecs9)

    LatVecs10 = gen_lattice_orthorhombic_P( 10.0, 3.0, 5.0 )
    println("orthorhombic")
    print_matrix(LatVecs10)

    LatVecs11 = gen_lattice_monoclinic_P( 10.0, 3.0, 5.0, 45.0 )
    println("monoclinic")
    print_matrix(LatVecs11)

    LatVecs12 = gen_lattice_triclinic( 10.0, 3.0, 5.0, 46.0, 45.0, 90.0 )
    println("triclinic")
    print_matrix(LatVecs12)

    return nothing
end

function test_atoms()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    atoms.LatVecs = gen_lattice_cubic(10.0)
    println(atoms)

    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """
    )
    atoms.LatVecs = gen_lattice_cubic(10.0)
    println(atoms)

    return nothing
end

@test test_gen_lattice() == nothing
@test test_atoms() == nothing

println("")
println("All test successfully run")
println("")

