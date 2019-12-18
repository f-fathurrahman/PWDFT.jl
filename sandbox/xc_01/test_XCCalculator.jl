using PWDFT

include("XCCalculator.jl")

function calc_xc( xc_calc::XCCalculator, Rhoe )
    println("Using XCCalculator")
    return
end

function calc_xc( xc_calc::LibxcXCCalculator, Rhoe )
    println("Using LibxcCalculator")
end

struct MyHamiltonian
    xc_calc::AbstractXCCalculator
end

function main()

    Ham1 = MyHamiltonian( LibxcXCCalculator() )
    println(Ham1)
    calc_xc( Ham1.xc_calc, 1.0 )

    Ham1 = MyHamiltonian( XCCalculator() )
    println(Ham1)
    calc_xc( Ham1.xc_calc, 1.0 )

end

main()