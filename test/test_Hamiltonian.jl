# H atom in a cube with pseudopotential
function test_Hamiltonian_v1()

    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        [joinpath(DIR_PSP, "pade_gth", "H-q1.gth")],
        15.0
    )

    return nothing
end

# H atom in a cube, all electron (without pseudopotential)
#function test_Hamiltonian_v2()
#    Ham = Hamiltonian(    
#        Atoms( xyz_string="""
#        1
#
#        H  0.0  0.0  0.0
#        """, LatVecs=gen_lattice_sc(16.0) ),
#        15.0
#    )
#
#    return nothing
#end


function test_Hamiltonian_v3()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "GBRV_LDA", "si_lda_v1.uspp.F.UPF")]
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, dual=dual, meshk=[3,3,3] )
    return nothing
end


function test_Hamiltonian_v4()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "PSLIB_US_PAW_LDA", "Si.pz-n-kjpaw_psl.1.0.0.UPF")]
    # FIXME: PSLIB_US_PAW_LDA is not included in the repo
    ecutwfc = 20.0 # or 40 Ry
    ecutrho = 100.0 # or 200 Ry
    dual = ecutrho/ecutwfc
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, dual=dual, meshk=[3,3,3] )
    return nothing
end


@testset "create Hamiltonian" begin
    @test test_Hamiltonian_v1() == nothing
    @test test_Hamiltonian_v3() == nothing
    @test test_Hamiltonian_v4() == nothing
end
