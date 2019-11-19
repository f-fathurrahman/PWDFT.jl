using PWDFT
using Libxc
using Printf

include("XC_x_slater.jl")
include("XC_c_pw.jl")
include("XC_c_vwn.jl")

include("XC_x_slater_spin.jl")
include("XC_c_pw_spin.jl")
include("XC_c_vwn_spin.jl")

include("XC_x_pbe.jl")
include("XC_c_pbe.jl")
include("XC_c_pbe_v2.jl")

include("XC_c_pbe_spin.jl")

include("Libxc_LDA.jl")
include("Libxc_GGA.jl")


function test_C_PBE()
    Rhoe = 1.2
    gRhoe2 = 0.1
    res1, _ = XC_c_pbe(Rhoe, gRhoe2)

    res4, _ = XC_c_pbe_v2(Rhoe, gRhoe2)

    res2 = calc_epsxc_GGA( [Rhoe], [gRhoe2], Libxc.GGA_C_PBE )
    res3 = calc_epsxc_LDA( [Rhoe], Libxc.LDA_C_PW )

    @printf("QE    C PBE v1 = %18.10f\n", res1)
    @printf("QE    C PBE v2 = %18.10f\n", res4)
    @printf("Libxc C PBE    = %18.10f\n", res2[1] - res3[1])

    @printf("QE    C PW = %18.10f\n", XC_c_pw(Rhoe)[1])
    @printf("Libxc C PW = %18.10f\n", res3[1])
end
test_C_PBE()

function test_X_PBE()
    Rhoe = 1.0
    gRhoe2 = 1.1
    res1, _ = XC_x_pbe(Rhoe, gRhoe2)

    res2 = calc_epsxc_GGA( [Rhoe], [gRhoe2], Libxc.GGA_X_PBE )
    res3 = calc_epsxc_LDA( [Rhoe], Libxc.LDA_X )

    @printf("QE    X PBE = %18.10f\n", res1)
    @printf("Libxc X PBE = %18.10f\n", res2[1] - res3[1])

    @printf("QE    X slater = %18.10f\n", XC_x_slater(Rhoe)[1])
    @printf("Libxc X slater = %18.10f\n", res3[1])
end
test_X_PBE()


function test_main()
    println( XC_x_slater(1.2) )
    println( XC_c_pw(1.2) )

    println( XC_c_vwn(1.2) )

    println( XC_x_slater_spin(1.2, 0.1) )
    println( XC_c_pw_spin(1.2, 0.1) )    
    println( XC_c_vwn_spin(1.2, 0.1) )

    println( XC_x_pbe(1.2, 1.1) )
    println( XC_c_pbe(1.2, 1.1) )

    println( XC_c_pbe_spin(1.2, 0.1, 1.1) )
end
#test_main()