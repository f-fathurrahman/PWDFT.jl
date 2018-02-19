using PWDFT
using PWDFT.PW

function test_main()
    LatVecs = 8.0*diagm(ones(3))
    pw = PWGrid(30.0, LatVecs, verbose=true)
    actual = prod(pw.Ns)/pw.gvecw.Ngwx
    theor = 1/(4*pi*0.25^3/3)
    @printf("Actual, theor: %10.5f %10.5f\n", actual, theor)
end

@time test_main()
@time test_main()
