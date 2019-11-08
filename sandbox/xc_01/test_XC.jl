using Printf

include("XC_x_slater.jl")
include("XC_c_pw.jl")
include("XC_c_vwn.jl")

include("XC_c_vwn_spin.jl")

include("XC_x_pbe.jl")
include("XC_c_pbe.jl")

function test_main()
    println( XC_x_slater(1.2) )
    println( XC_c_pw(1.2) )

    println( XC_c_vwn(1.2) )

    println( XC_c_vwn_spin(1.2, 0.1) )

    println( XC_x_pbe(1.2, 1.1) )
    println( XC_c_pbe(1.2, 1.1) )
end

test_main()