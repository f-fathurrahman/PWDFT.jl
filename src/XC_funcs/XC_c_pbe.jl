# PBE correlation (without LDA part)
# iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
function XC_c_pbe( rho, grho )

    ga = 0.0310906908696548950
    be = 0.06672455060314922
    third = 1.0/3.0
    pi34 = 0.6203504908994
    xkf = 1.919158292677513
    xks = 1.128379167095513

    rs = pi34/rho^third
    ec, vc = XC_c_pw( rho )
    
    kf = xkf/rs
    ks = xks * sqrt(kf)
    t = sqrt(grho) / (2.0 * ks * rho)
    
    expe = exp(-ec/ga)
    af = be / ga * (1.0 / (expe - 1.0) )
    
    bf = expe * (vc - ec)
  
    y = af * t * t
  
    xy = (1.0 + y) / (1.0 + y + y * y)
  
    qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
  
    s1 = 1.0 + be / ga * t * t * xy
    h0 = ga * log(s1)
  
    dh0 = be * t * t / s1 * ( -7.0 / 3.0 * xy - qy * (af * bf / be - 7.0 / 3.0) )
  
    ddh0 = be/(2.0 * ks * ks * rho) * (xy - qy) / s1
  
    sc = rho * h0
    v1c = h0 + dh0
    v2c = ddh0

    return sc, v1c, v2c

end