using Printf
using PWDFT

include("calc_ylm_angles.jl")
include("Ylm_real_v2.jl")

function test_main()

    g = [0.0, 1.0, 2.0]
    Ng = 1
    lmax = 3

    calc_ylm_angles(g)

    lm = 0
    for l = 0:lmax
        
        @printf("\nl = %d\n", l)
        
        m = 0
        @printf("m = %d\n", m)
        lm = lm + 1
        ylm = Ylm_real(l, m, g)
        @printf("%18.10f\n", ylm)

        for m = 1:l
            
            @printf("m = %d\n", m)
            lm = lm + 1
            ylm = Ylm_real(l, m, g)
            @printf("%18.10f\n", ylm)
            
            @printf("m = %d\n", -m)
            lm = lm + 1
            ylm = Ylm_real(l, -m, g)
            @printf("%18.10f\n", ylm)

        end
    end

    println("lm = ", lm)

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


#test_second()


