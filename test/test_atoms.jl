function test_atoms()

    atoms = init_atoms_xyz(
        joinpath( DIR_PWDFT, "structures/CuSO4.xyz"),
        verbose=true
    )
    atoms.LatVecs = gen_lattice_cubic(10.0)

    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """,
        verbose=true
    )
    atoms.LatVecs = gen_lattice_cubic(10.0)

    return nothing
end

@testset "init_atoms function" begin
    @test test_atoms() == nothing
end