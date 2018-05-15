using PWDFT
using Base.Test

function test_gen_lattice()
    return true
end

@test gen_lattice_cubic(16.0) == 16.0*diagm(ones(3))

println("All test successfully run")



