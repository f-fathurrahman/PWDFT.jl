using Printf
using SpecialFunctions
using QuadGK

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()

    #psp = PsPot_GTH(joinpath(DIR_PSP, "H-q1.gth"))
    psp = PsPot_GTH(joinpath(DIR_PSP, "Si-q4.gth"))
    println(psp)

    Zval = psp.zval

    Gm = 10.0

    myfunc(r) = r^2*besselj0(Gm*r)*( PWDFT.eval_Vloc_R( psp, r ) + Zval*erf(r)/r )
    #myfunc(r) = r^2*( PWDFT.eval_Vloc_R( psp, r ) + Zval*erf(r)/r )

    
    for rmax in [1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0, 1000.0]
        res = quadgk( myfunc, eps(), rmax )
        @printf("rmax = %18.10f, integ = %18.10f, err = %18.10e\n", rmax, res[1], res[2])
    end

end

main()
