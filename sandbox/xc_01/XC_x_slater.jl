# Adapted from PWSCF
function XC_x_slater( rs )
    f = -0.687247939924714
    alpha = 2.0/3.0

    ex = f * alpha / rs
    vx = 4.0 / 3.0 * f * alpha / rs
    return ex, vx
end