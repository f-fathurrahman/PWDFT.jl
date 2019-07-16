atoms = Atoms(xyz_string_frac=
    """
    2

    Si  0.0  0.0  0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true)

@testset "Kpoints.jl" begin
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions

    kpoints = KPoints( atoms, [3,3,3], [1,1,1] )
    println(kpoints)
    # QE and VASP output 6 ir_kpoints
    # while here only 4 ir_kpoints
    # From symmetry perspective, only 4 kpoints are irreducible
    # Get the QE and Vasp result by setting atom position following(break symmetry):
    # Si  0.0  0.0  0.0
    # Si  0.26  0.26  0.26

    # test KPoints constructor dispatch
    kpoints = KPoints( atoms, (3,3,3), (1,1,1) )
    kpoints = KPoints( atoms, [3,3,3], (1,1,1) )
    kpoints = KPoints( atoms, (3,3,3), [1,1,1] )
end
