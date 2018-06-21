using Printf
using PWDFT

function test_main()

    pw = PWGrid( 5.0, gen_lattice_sc(3.0) )
    G = pw.gvec.G

    println(pw)

    l = 1
    m = -1

    Ng = size(G)[2]
    println("Ng = ", Ng)

    for ig = 1:Ng
        g = G[:,ig]
        ylm = Ylm_real(l, m, g)
        @printf("%18.10f %18.10f %18.10f %18.10f\n", g[1], g[2], g[3], ylm)
    end

end

function test_second()
    l = 1
    m = 1
 
    g = [1.0, 0.0, 0.0]
    ylm = Ylm_real(l,m,g)
    @printf("Ylm = %18.10f\n", ylm)

    g = [0.0, 1.0, 0.0]
    ylm = Ylm_real(l,m,g)
    @printf("Ylm = %18.10f\n", ylm)

    g = [0.0, 0.0, 1.0]
    ylm = Ylm_real(l,m,g)
    @printf("Ylm = %18.10f\n", ylm)

end


test_main()
test_second()


