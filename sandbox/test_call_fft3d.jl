using LinearAlgebra
using Random
using Printf
using PWDFT

function test_multicolumn( ecutwfc_Ry::Float64, Nstates::Int64 )

    pw = PWGrid( 0.5*ecutwfc_Ry, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints, Nstates )

    println("\n\nMulti column version")
    println("Data size: ", Ns)
    println("Nstates  : ", Nstates)

    println("\nC version")
    @time out1c = c_R_to_G( Ns, in1 )
    @time in2c = c_G_to_R( Ns, out1c )

    println("\nUsing plan")
    @time out2 = R_to_G( pw, in1 )
    @time in3 = G_to_R( pw, out2 )

    println("\nNaive version")
    @time out1 = R_to_G( Ns, in1 )
    @time in2 = G_to_R( Ns, out1 )

    println("\nCheck")
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))

    println("diff out2 = ", sum( abs.(out1c - out2) ))
    println("diff out2 = ", sum( abs.(in2c - in3) ))

end



function test_singlecolumn( ecutwfc_Ry::Float64 )

    pw = PWGrid( 0.5*ecutwfc_Ry, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    in1 = rand(ComplexF64, Npoints)

    println("\n\nSingle column version")
    println("Data size:", Ns)

    println("\nC version")
    @time out1c = c_R_to_G( Ns, in1 )
    @time in2c = c_G_to_R( Ns, out1c )

    println("\nUsing plan")
    @time out2 = R_to_G( pw, in1 )
    @time in3 = G_to_R( pw, out2 )

    println("\nNaive version")
    @time out1 = R_to_G( Ns, in1 )
    @time in2 = G_to_R( Ns, out1 )

    println("\nCheck")
    println("diff out1 = ", sum( abs.(out1c - out1) ))
    println("diff out1 = ", sum( abs.(in2c - in2) ))

    println("diff out2 = ", sum( abs.(out1c - out2) ))
    println("diff out2 = ", sum( abs.(in2c - in3) ))

end

test_singlecolumn( 30.0 )
test_singlecolumn( 40.0 )
test_singlecolumn( 50.0 )

test_multicolumn( 30.0, 10 )
test_multicolumn( 40.0, 10 )
test_multicolumn( 50.0, 10 )
