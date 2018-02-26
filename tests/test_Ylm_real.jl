using PWDFT

function test_main()

  pw = PWGrid( 5.0, 3.0*diagm(ones(3)) )
  G = pw.gvec.G

  println(pw)

  l = 1
  m = 0

  Ng = size(G)[2]
  println("Ng = ", Ng)
  exit()

  for ig = 1:Ng
      g = G[:,ig]
      ylm = Ylm_real(l, m, g)
      @printf("%18.10f %18.10f %18.10f %18.10f\n", g[1], g[2], g[3], ylm)
  end

end

test_main()
