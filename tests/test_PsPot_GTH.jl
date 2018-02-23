using PWDFT

function test_main()
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q10.gth")
    println(psp)
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q18.gth")
    println(psp)
end

test_main()
