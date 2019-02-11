using Printf
using PWDFT
using BenchmarkTools

function test_main()
    energies = Energies()
    println(energies)

    println("test sum: ", sum(energies))
end

function test_bench()
    energies = Energies()
    @btime Etot = sum($(energies))
end

test_main()
test_bench()


