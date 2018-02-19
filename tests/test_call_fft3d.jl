using PWDFT

function test_multicolumn(Ns)
    Nstates = 4
    in1 = rand(Complex128,prod(Ns),Nstates)

    @time out1c = c_R_to_G(Ns,in1)
    @time in2c = c_G_to_R(Ns,out1c)

    @time out1 = R_to_G(Ns,in1)
    @time in2 = G_to_R(Ns,out1)

    println("\nMulticolumn")
    println("Data size:", Ns)
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))
end

function test_singlecolumn(Ns)
    in1 = rand(Complex128,prod(Ns))

    @time out1c = c_R_to_G(Ns,in1)
    @time in2c = c_G_to_R(Ns,out1c)

    @time out1 = R_to_G(Ns,in1)
    @time in2 = G_to_R(Ns,out1)

    println("\nSingle column")
    println("Data size:", Ns)
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))
end

test_singlecolumn([5,5,5])
test_multicolumn([5,5,5])

test_singlecolumn([45,45,45])
test_multicolumn([45,45,45])

