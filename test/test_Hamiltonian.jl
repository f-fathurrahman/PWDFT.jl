# H atom in a cube with pseudopotential
function test_Hamiltonian_v1()
    
    dir_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        [joinpath(dir_PWDFT, "pseudopotentials/pade_gth/H-q1.gth")],
        15.0
    )
    println(Ham)
    return nothing
end

# H atom in a cube, all electron (without pseudopotential)
function test_Hamiltonian_v2()
    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        15.0
    )
    println(Ham)
    return nothing
end