using PWDFT

function test_main()
    psp = PsPot_HGH("H", "../LDA_HGH/H.1.hgh")
    println(psp)

    psp = PsPot_HGH("Ge", "../LDA_HGH/Ge.4.hgh")
    println(psp)
end

test_main()
