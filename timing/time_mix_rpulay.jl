using PWDFT
using Printf
using BenchmarkTools
using Random

Random.seed!(1234)

function do_mix_rpulay(Rhoe, Rhoe_new, XX, FF, x_old, f_old)
    Npoints = size(Rhoe)
    MIXDIM = size(XX,2)
    betamix = 0.1
    iter = 6
    Rhoe = mix_rpulay!(Rhoe, Rhoe_new, betamix, XX, FF, iter, MIXDIM, x_old, f_old)
    return
end

function time_mix_rpulay()

    @printf("\n")
    @printf("-----------------\n")
    @printf("Timing mix_rpulay\n")
    @printf("-----------------\n")
    @printf("\n")

    NN = (20,30,40,50,60,70,80)
    MIXDIM = 5
    for N in NN 
        Npoints = N^3
        # create mock data
        Rhoe = rand(Npoints)
        Rhoe_new = rand(Npoints)
        XX = rand(Npoints,MIXDIM)
        FF = rand(Npoints,MIXDIM)
        x_old = rand(Npoints)
        f_old = rand(Npoints)
        @printf("Npoints = %8d: ", Npoints)
        @time do_mix_rpulay(Rhoe, Rhoe_new, XX, FF, x_old, f_old)
        @time do_mix_rpulay(Rhoe, Rhoe_new, XX, FF, x_old, f_old)
        println()
    end
end

time_mix_rpulay()

