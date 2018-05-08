using PWDFT

function test_main()
    LatVecs = 16.0*diagm(ones(3))
    ecutwfc_Ry = 30.0
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)
end

@time test_main()
@time test_main()
