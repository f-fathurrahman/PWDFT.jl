using PWDFT

function test_multicolumn(Ns)
    Nstates = 4
    in1 = rand(Complex128,prod(Ns),Nstates)

    println("\n\nMulticolumn version")
    println("Data size:", Ns)

    println("\nC version")
    @time out1c = c_R_to_G(Ns,in1)
    @time in2c = c_G_to_R(Ns,out1c)

    println("\nNaive version")
    @time out1 = R_to_G(Ns,in1)
    @time in2 = G_to_R(Ns,out1)

    println("\nCheck")
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))
end

function test_singlecolumn(Ns)
    Npoints = prod(Ns)
    in1 = rand(Complex128, Npoints)

    println("\n\nSingle column version")
    println("Data size:", Ns)

    println("\nC version")
    @time out1c = c_R_to_G(Ns,in1)
    @time in2c = c_G_to_R(Ns,out1c)

    println("\nNaive version")
    @time out1 = R_to_G(Ns,in1)
    @time in2 = G_to_R(Ns,out1)

    println("\nCheck")
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))

    println("\nTime for creating plan")
    @time planfw = init_plan_forward(Ns)
    println("\nTime for executing plan")
    @time out1p = reshape( planfw*reshape(in1, Ns[1], Ns[2], Ns[3]), Npoints )

    println("diff with planfw = ", sum( abs.(out1 - out1p) ))

end

function R_to_G( planfw, Ns, fR::Array{Complex128,1} )
    return
end

function init_plan_forward( Ns::Array{Int64} )
    return plan_fft( zeros(Ns[1],Ns[2],Ns[3]) )
end



test_singlecolumn([5,5,5])
test_multicolumn([5,5,5])

test_singlecolumn([45,45,45])
test_multicolumn([45,45,45])

test_singlecolumn([60,60,60])
test_multicolumn([60,60,60])
