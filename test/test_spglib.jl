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
    num_ir, _, _ = spg_get_ir_reciprocal_mesh(atoms, mesh, is_shift)
    @test num_ir == 8
end
