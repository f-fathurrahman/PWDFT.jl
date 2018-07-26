using Printf

include("../src/Energies.jl")

function test_main()
    energies = Energies()
    println(energies)

    println("test sum: ", sum(energies))
end

using BenchmarkTools

function test_bench()
    energies = Energies()
    @btime Etot = sum($(energies))
end

test_main()
test_bench()


