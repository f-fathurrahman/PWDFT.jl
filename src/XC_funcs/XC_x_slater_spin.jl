#
# Slater exchange with alpha=2/3, spin-polarized case
#
function XC_x_slater_spin( rho, zeta )

    f = -1.10783814957303361 # f = -9/8*(3/pi)^(1/3)
    alpha = 2.0/3.0
    third = 1.0/3.0
    p43 = 4.0/3.0

    rho13 = ( (1.0 + zeta) * rho)^third
    exup = f * alpha * rho13
    vxup = p43 * f * alpha * rho13
    
    rho13 = ( (1.0 - zeta) * rho)^third
    exdw = f * alpha * rho13
    vxdw = p43 * f * alpha * rho13
    ex = 0.5 * ( (1.0 + zeta) * exup + (1.0 - zeta) * exdw)

    return ex, vxup, vxdw
end