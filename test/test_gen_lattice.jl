function test_gen_lattice()

    LatVecs1 = gen_lattice_cubic(10.0)

    LatVecs2 = gen_lattice_fcc(10.0)

    LatVecs3 = gen_lattice_bcc(10.0)
    LatVecs4 = gen_lattice_bcc_v2(10.0)
    
    @assert det(LatVecs3) == det(LatVecs4)

    LatVecs5 = gen_lattice_hexagonal( 10.0, 15.0 )

    LatVecs6 = gen_lattice_trigonal( 10.0, 90.0 )

    LatVecs7 = gen_lattice_trigonal_v2( 10.0, 90.0 )

    LatVecs8 = gen_lattice_tetragonal_P( 10.0, 5.0 )

    LatVecs9 = gen_lattice_tetragonal_I( 10.0, 5.0 )

    LatVecs10 = gen_lattice_orthorhombic( 10.0, 3.0, 5.0 )

    LatVecs11 = gen_lattice_monoclinic( 10.0, 3.0, 5.0, 45.0 )

    LatVecs12 = gen_lattice_triclinic( 10.0, 3.0, 5.0, 46.0, 45.0, 90.0 )

    return nothing
end


@testset "gen_lattice" begin
    @test gen_lattice_cubic(16.0) == 16.0*Matrix(Diagonal(ones(3)))
    @test test_gen_lattice() == nothing
end