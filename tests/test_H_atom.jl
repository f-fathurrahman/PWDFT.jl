using PWDFT
using PWDFT.PW


include("init_V_coulomb_G.jl")
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
    #
    V = init_V_coulomb_G( pw, strf, [1.0] )
    println("sum V = ", sum(V))
end

test_main()
