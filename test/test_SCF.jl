# H atom in a cube with pseudopotential
function test_SCF_v1()
    
    dir_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        [joinpath(dir_PWDFT, "pseudopotentials/pade_gth/H-q1.gth")],
        15.0
    )
    
    KS_solve_SCF!( Ham, NiterMax=1 )
    KS_solve_SCF!( Ham, betamix=0.5, NiterMax=1  )
    KS_solve_SCF!( Ham, mix_method="pulay", NiterMax=6 )
    KS_solve_SCF!( Ham, mix_method="ppulay", NiterMax=6 )
    KS_solve_SCF!( Ham, mix_method="rpulay", NiterMax=6 )
    KS_solve_SCF!( Ham, mix_method="anderson", NiterMax=6 )

    return nothing
end

# H atom in a cube, all electron (without pseudopotential)
function test_SCF_v2()
    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        15.0
    )
    #KS_solve_SCF!( Ham )
    KS_solve_SCF!( Ham, mix_method="pulay", NiterMax=6 )

    return nothing
end