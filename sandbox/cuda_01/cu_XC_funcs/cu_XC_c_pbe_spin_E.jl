# epsxc only version
function cu_XC_c_pbe_spin_E( rho, zeta, grho )

    ga = 0.031091
    be = 0.06672455060314922
    third = 1.0/3.0
    pi34 = 0.6203504908994
  
    xkf = 1.919158292677513
    xks = 1.128379167095513

    rs = pi34/CUDAnative.pow(rho, third)
    
    ec = cu_XC_c_pw_spin_E( rho, zeta )

    kf = xkf / rs
    ks = xks * CUDAnative.sqrt(kf)
    fz = 0.5 * ( CUDAnative.pow( 1.0 + zeta, 2.0/3.0 ) + CUDAnative.pow( 1.0 - zeta, 2.0/3.0) )
    
    fz2 = fz * fz
    fz3 = fz2 * fz
    
    t = CUDAnative.sqrt(grho) / (2.0 * fz * ks * rho)
    expe = CUDAnative.exp( -ec / (fz3 * ga) )
    af = be / ga * (1.0 / (expe-1.0) )
    
    y = af * t * t
    
    xy = (1.0 + y) / (1.0 + y + y * y)
    qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
    
    s1 = 1.0 + be / ga * t * t * xy
    h0 = fz3 * ga * CUDAnative.log(s1)
  
    sc = rho * h0
  
    return sc

end