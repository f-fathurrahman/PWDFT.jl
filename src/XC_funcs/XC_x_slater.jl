function XC_x_slater( Rhoe::Float64 )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    f = -0.687247939924714
    alpha = 2.0/3.0

    ex = f * alpha / rs
    vx = 4.0 / 3.0 * f * alpha / rs
    return ex, vx
end
