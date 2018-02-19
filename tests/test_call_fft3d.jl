using PWDFT

function test_main()
  Ns = [30,30,30]
  in1 = rand(Complex128,prod(Ns))

  @time out1c = c_R_to_G(Ns,in1)
  @time in2c = c_G_to_R(Ns,out1c)

  @time out1 = R_to_G(Ns,in1)
  @time in2 = G_to_R(Ns,out1)

  println("diff out1 = ", sum( abs.(out1c - out1) ))
  println("diff out1 = ", sum( abs.(in2c - in2) ))
end

test_main()

