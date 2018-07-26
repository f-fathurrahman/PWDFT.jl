using Printf

include("../src/Energies.jl")

function test_main()
    energies = Energies()
    println(energies)

    println("test sum: ", sum(energies))
end

test_main()
