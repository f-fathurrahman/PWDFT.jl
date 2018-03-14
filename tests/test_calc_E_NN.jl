using PWDFT


include("../src/calc_E_NN_simple.jl")
include("../src/calc_E_NN_new.jl")


"""
ewald = 3.58126919 Ry = 1.790634595 Ha
"""
function test_N2()
    atoms = init_atoms_xyz("N2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [5.0] # XXX hardwired

    E_NN = calc_E_NN_simple( LatVecs, atoms, Zvals, verbose=true )
    E_NN_v2 = calc_E_NN( PWGrid(25.0,LatVecs), atoms, Zvals )
    E_NN_v3 = calc_E_NN_new( LatVecs, atoms, Zvals ) 

    E_NN_ref = 1.790634595
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)
    d3 = abs(E_NN_v3 - E_NN_ref)

    @printf("\n version simple, now, and new\n\n")
    @printf("%18.10f %18.10f %18.10f\n", E_NN, E_NN_v2, E_NN_v3)
    @printf("%18.7e %18.7e %18.7e\n", d1, d2, d3)

end

"""
ewald contribution = 0.62633998 Ry = 0.31316999 Ha
"""
function test_H2()
    atoms = init_atoms_xyz("H2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN_simple( LatVecs, atoms, Zvals, verbose=true )
    E_NN_v2 = calc_E_NN( PWGrid(25.0,LatVecs), atoms, Zvals )
    E_NN_v3 = calc_E_NN_new( LatVecs, atoms, Zvals )    

    E_NN_ref = 0.31316999 
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)
    d3 = abs(E_NN_v3 - E_NN_ref)

    @printf("\n version simple, now, and new\n\n")
    @printf("%18.10f %18.10f %18.10f\n", E_NN, E_NN_v2, E_NN_v3)
    @printf("%18.7e %18.7e %18.7e\n", d1, d2, d3)

end

"""
PWSCF result = -0.04393722 Ry = -0.02196861 Ha
"""
function test_LiH()
    atoms = init_atoms_xyz("LiH.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0,1.0] # XXX hardwired

    E_NN = calc_E_NN_simple( LatVecs, atoms, Zvals, verbose=true )
    E_NN_v2 = calc_E_NN( PWGrid(25.0,LatVecs), atoms, Zvals )
    E_NN_v3 = calc_E_NN_new( LatVecs, atoms, Zvals )

    E_NN_ref = -0.02196861
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)
    d3 = abs(E_NN_v3 - E_NN_ref)

    @printf("\n version simple, now, and new\n\n")
    @printf("%18.10f %18.10f %18.10f\n", E_NN, E_NN_v2, E_NN_v3)
    @printf("%18.7e %18.7e %18.7e\n", d1, d2, d3)

end

"""
PWSCF result = -0.17733109 Ry = -0.088665545 Ha
"""
function test_H()
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN_simple( LatVecs, atoms, Zvals, verbose=true )
    E_NN_v2 = calc_E_NN( PWGrid(25.0,LatVecs), atoms, Zvals )
    E_NN_v3 = calc_E_NN_new( LatVecs, atoms, Zvals )

    E_NN_ref = -0.088665545
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)
    d3 = abs(E_NN_v3 - E_NN_ref)

    @printf("\n version simple, now, and new\n\n")
    @printf("%18.10f %18.10f %18.10f\n", E_NN, E_NN_v2, E_NN_v3)
    @printf("%18.7e %18.7e %18.7e\n", d1, d2, d3)

end

#test_H()
#test_H2()
#test_LiH()
test_N2()

