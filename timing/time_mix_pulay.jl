using PWDFT

function time_mix_pulay()

    Npoints = 500
    Rhoe = rand(Npoints)
    Rhoe_new = rand(Npoints)

    MIXDIM = 5
    XX = rand(Npoints,MIXDIM)
    FF = rand(Npoints,MIXDIM)

    x_old = rand(Npoints)
    f_old = rand(Npoints)

    betamix = 0.1

    iter = 6

    # create mock data
    Rhoe = mix_pulay!(Rhoe, Rhoe_new, betamix, XX, FF, iter, MIXDIM, x_old, f_old)

end

time_mix_pulay()

