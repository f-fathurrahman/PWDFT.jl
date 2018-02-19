using PWDFT
using PWDFT.PW

function test_main()
    LatVecs = 16.0*diagm(ones(3))
    pw = PWGrid(30.0, LatVecs)
    println(pw)
end

@time test_main()
@time test_main()
