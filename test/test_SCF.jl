# H atom in a cube with pseudopotential
function test_SCF_v1()
    Ham = Hamiltonian(    
        Atoms( xyz_string="""
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0) ),
        ["../pseudopotentials/pade_gth/H-q1.gth"],
        15.0
    )
    
    KS_solve_SCF!( Ham )
    KS_solve_SCF!( Ham, betamix=0.5 )
    KS_solve_SCF!( Ham, mix_method="pulay" )
    KS_solve_SCF!( Ham, mix_method="ppulay" )
    KS_solve_SCF!( Ham, mix_method="rpulay" )
    KS_solve_SCF!( Ham, mix_method="anderson" )

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
    KS_solve_SCF!( Ham, mix_method="pulay" )

    return nothing
end