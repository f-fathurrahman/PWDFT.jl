atoms = Atoms(xyz_string_frac=
    """
    2

    Si  0.0  0.0  0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

@testset "spglib.jl" begin
    n, rots, trans = spg_get_symmetry(atoms)
    @test n == 48
    @test size(rots) == (3, 3, 48)
    @test size(trans) == (3, 48)


    mesh = [4, 4, 4]
    is_shift = [0, 0, 0]
    num_ir, kgrid, mapping = spg_get_ir_reciprocal_mesh(atoms, mesh, is_shift)
    @test num_ir == 8
    @test kgrid == [0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1 0 1 2 -1; 0 0 0 0 1 1 1 1 2 2 2 2 -1 -1 -1 -1 0 0 0 0 1 1 1 1 2 2 2 2 -1 -1 -1 -1 0 0 0 0 1 1 1 1 2 2 2 2 -1 -1 -1 -1 0 0 0 0 1 1 1 1 2 2 2 2 -1 -1 -1 -1; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    @test mapping == [0, 1, 2, 1, 1, 5, 6, 7, 2, 6, 10, 6, 1, 7, 6, 5, 1, 5, 6, 7, 5, 1, 7, 6, 6, 7, 6, 27, 7, 6, 27, 6, 2, 6, 10, 6, 6, 7, 6, 27, 10, 6, 2, 6, 6, 27, 6, 7, 1, 7, 6, 5, 7, 6, 27, 6, 6, 27, 6, 7, 5, 6, 7, 1]
end
