using PWDFT

include("calc_E_NN.jl")
include("calc_strfact.jl")

function test_main()
    const LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0*0.5
    pw = PWGrid( ecutwfc_Ry, LatVecs )
    println(pw)
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)
    #
    strf = calc_strfact( atoms, pw )

    E_NN = calc_E_NN( pw, strf, atoms.positions, atoms.Nspecies, atoms.atm2species, [1.0])
end

test_main()
