# PBE exchange (without Slater exchange):
# iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
function XC_x_pbe( rho, grho )

    # From xc_gga_drivers.f90
    SMALL1 = 1e-6 # rho_threshold
    SMALL2 = 1e-10 # grho_threshold
    if rho <= SMALL1
        return 0.0, 0.0, 0.0
    end
    if grho <= SMALL2
        return 0.0, 0.0, 0.0
    end

    third = 1.0/3.0
    k = 0.804
    mu = 0.2195149727645171
    c1 = 0.75/pi
    c2 = 3.093667726280136
    c5 = 4.0*third

    agrho = sqrt(grho)
    kf = c2 * rho^third
    dsg = 0.5/kf
    s1 = agrho*dsg / rho
    s2 = s1 * s1
    ds = -c5 * s1
    
    # Energy
    f1 = s2 * mu/k
    f2 = 1.0 + f1
    f3 = k / f2
    fx = k - f3
    exunif = -c1 * kf
 
    sx_s = exunif * fx
    
    dxunif = exunif * third
    ds = -c5 * s1
    
    dfx1 = f2 * f2
    dfx = 2.0 * mu * s1 / dfx1
    
    v1x = sx_s + dxunif * fx + exunif * dfx * ds
    v2x = exunif * dfx * dsg / agrho
    sx  = sx_s * rho

    return sx, v1x, v2x
end