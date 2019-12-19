#
# PBE correlation (without LDA part) - spin-polarized
# iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
#
function XC_c_pbe_spin( rho, zeta, grho )
  
    ga = 0.031091
    be = 0.06672455060314922
    third = 1.0/3.0
    pi34 = 0.6203504908994
  
    xkf = 1.919158292677513
    xks = 1.128379167095513

    rs = pi34/rho^third
    
    ec, vcup, vcdw = XC_c_pw_spin( rho, zeta )

    kf = xkf / rs
    ks = xks * sqrt(kf)
    fz = 0.5 * ( (1.0 + zeta)^(2.0 / 3.0) + (1.0 - zeta)^(2.0 / 3.0) )
    
    fz2 = fz * fz
    fz3 = fz2 * fz
    fz4 = fz3 * fz
    dfz = ( (1.0 + zeta)^(-1.0/3.0) - (1.0 - zeta)^( -1.0 / 3.0) ) / 3.0
    
    t = sqrt(grho) / (2.0 * fz * ks * rho)
    expe = exp( - ec / (fz3 * ga) )
    af = be / ga * (1.0 / (expe-1.0) )
    
    bfup = expe * (vcup - ec) / fz3
    bfdw = expe * (vcdw - ec) / fz3
    
    y = af * t * t
    
    xy = (1.0 + y) / (1.0 + y + y * y)
    qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
    
    s1 = 1.0 + be / ga * t * t * xy
    h0 = fz3 * ga * log(s1)
    
    dh0up = be * t * t * fz3 / s1 * ( -7.0 / 3.0 * xy - qy * (af * bfup / be - 7.0 / 3.0) )
    
    dh0dw = be * t * t * fz3 / s1 * ( -7.0 / 3.0 * xy - qy * (af * bfdw / be - 7.0 / 3.0) )
    
    dh0zup = (3.0 * h0 / fz - be * t * t * fz2 / s1 * (2.0 * xy - 
             qy * (3.0 * af * expe * ec / fz3 / be + 2.0) ) ) * dfz * (1.0 - zeta)
    
    dh0zdw = -(3.0 * h0 / fz - be * t * t * fz2 / s1 * (2.0 * xy -
             qy * (3.0 * af * expe * ec / fz3 / be + 2.0) ) ) * dfz * (1.0 + zeta)

    ddh0 = be * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1
  
    sc = rho * h0
    v1cup = h0 + dh0up + dh0zup
    v1cdw = h0 + dh0dw + dh0zdw
    v2c = ddh0
  
    return sc, v1cup, v1cdw, v2c

end