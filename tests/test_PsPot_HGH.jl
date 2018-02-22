using PWDFT

function test_main()
    psp = PsPot_HGH("Ge", "../LDA_HGH/Ge.4.hgh")
    println(psp)
end

test_main()
