# energy (epsxc) only version
function XC_x_slater_spin_E( rho, zeta )

    f = -1.10783814957303361 # f = -9/8*(3/pi)^(1/3)
    alpha = 2.0/3.0
    third = 1.0/3.0
    p43 = 4.0/3.0

    rho13 = ( (1.0 + zeta) * rho)^third
    exup = f * alpha * rho13
    
    rho13 = ( (1.0 - zeta) * rho)^third
    exdw = f * alpha * rho13

    ex = 0.5 * ( (1.0 + zeta) * exup + (1.0 - zeta) * exdw)

    return ex
end