using PWDFT

"""
ewald = 3.58126919 Ry = 1.790634595 Ha
"""
function test_N2()
    atoms = init_atoms_xyz("N2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [5.0] # XXX hardwired

    E_NN = calc_E_NN( LatVecs, atoms, Zvals )

    E_NN_ref = 1.790634595
    diffE = abs(E_NN - E_NN_ref)

    @printf("E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
ewald contribution = 0.62633998 Ry = 0.31316999 Ha
"""
function test_H2()
    atoms = init_atoms_xyz("H2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN( LatVecs, atoms, Zvals )

    E_NN_ref = 0.31316999 
    diffE = abs(E_NN - E_NN_ref)

    @printf("E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
PWSCF result = -0.04393722 Ry = -0.02196861 Ha
"""
function test_LiH()
    atoms = init_atoms_xyz("LiH.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0,1.0] # XXX hardwired

    E_NN = calc_E_NN( LatVecs, atoms, Zvals )
    
    E_NN_ref = -0.02196861
    diffE = abs(E_NN - E_NN_ref)

    @printf("E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
PWSCF result = -0.17733109 Ry = -0.088665545 Ha
"""
function test_H()
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN( LatVecs, atoms, Zvals )

    E_NN_ref = -0.088665545
    diffE = abs(E_NN - E_NN_ref)

    @printf("E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

test_H()
test_H2()
test_LiH()
test_N2()

