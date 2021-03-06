using PWDFT
using Printf
using BenchmarkTools
using Random

Random.seed!(1234)

function do_mix_pulay(Npoints)

    # create mock data
    Rhoe = rand(Npoints)
    Rhoe_new = rand(Npoints)

    MIXDIM = 5
    XX = rand(Npoints,MIXDIM)
    FF = rand(Npoints,MIXDIM)

    x_old = rand(Npoints)
    f_old = rand(Npoints)

    betamix = 0.1

    iter = 6

    Rhoe = mix_pulay!(Rhoe, Rhoe_new, betamix, XX, FF, iter, MIXDIM, x_old, f_old)

end

function time_mix_pulay()

    @printf("\n")
    @printf("---------------------------------------------------------\n")
    @printf("Timing mix_pulay (includes time for allocating mock data)\n")
    @printf("---------------------------------------------------------\n")
    @printf("\n")

    NN = (20,30,40,50,60,70,80)
    for N in NN 
        Npoints = N^3
        @printf("Npoints = %8d: ", Npoints)
        @btime do_mix_pulay($Npoints)
    end
end

time_mix_pulay()

