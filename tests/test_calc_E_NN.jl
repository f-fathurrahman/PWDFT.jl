using PWDFT


include("../src/calc_E_NN_simple.jl")


function test_conv_E_nn_N2()
    atoms = init_atoms_xyz("N2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [5.0] # XXX hardwired

    for i = 1:10
        ecutwfc_Ry = 30.0 + (i-1)*10.0
        pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
        E_NN = calc_E_NN( pw, atoms, Zvals )
        E_NN_old = calc_E_NN_simple( pw, atoms, Zvals )
        @printf("%18.10f %18.10f %18.10f\n", ecutwfc_Ry, E_NN_old, E_NN)
    end

end

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

    E_NN_ref = 1.790634595
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)

    @printf("%18.10f %18.10f %18.10e %18.10e\n", E_NN, E_NN_v2, d1, d2)
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

    E_NN_ref = 0.31316999 
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)

    @printf("%18.10f %18.10f %18.10e %18.10e\n", E_NN, E_NN_v2, d1, d2)

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

    E_NN_ref = -0.02196861
    d1 = abs(E_NN - E_NN_ref)
    d2 = abs(E_NN_v2 - E_NN_ref)

    @printf("%18.10f %18.10f %18.10e %18.10e\n", E_NN, E_NN_v2, d1, d2)

end



function test_conv_E_nn_LiH()
    atoms = init_atoms_xyz("N2.xyz")
    println(atoms)

    LatVecs = 16.0*diagm( ones(3) )
    Zvals = [5.0] # XXX hardwired

    for i = 1:10
        ecutwfc_Ry = 30.0 + (i-1)*10.0
        pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
        E_NN = calc_E_NN( pw, atoms, Zvals )
        E_NN_old = OLD_calc_E_NN( pw, atoms, Zvals, verbose=true )
        @printf("%18.10f %18.10f %18.10f\n", ecutwfc_Ry, E_NN_old, E_NN)
    end

end

function test_main()
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 60.0
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)
    #
    strf = calc_strfact( atoms, pw )

    E_NN = OLD_calc_E_NN( pw, atoms.positions, atoms.Nspecies, atoms.atm2species, [1.0])
    @printf("E_NN = %18.10f\n", E_NN)
    
    E_NN = calc_E_NN( pw, atoms, [1.0] )
    @printf("E_NN = %18.10f\n", E_NN)
end

#test_main()
#test_conv_E_nn_N2()
test_N2()
#test_H2()
#test_LiH()

