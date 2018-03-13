using PWDFT


include("../src/OLD_calc_E_NN.jl")

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

    E_NN = OLD_calc_E_NN( pw, strf, atoms.positions, atoms.Nspecies, atoms.atm2species, [1.0])
    @printf("E_NN = %18.10f\n", E_NN)
    
    E_NN = calc_E_NN( pw, atoms, [1.0] )
    @printf("E_NN = %18.10f\n", E_NN)
end

test_main()
