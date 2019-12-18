# Adapted from PWSCF
# I found slight difference with Libxc for PBE correlation.

# PBE correlation (without LDA part)
# iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).

function cu_XC_c_pbe( rho, grho )

    ga = 0.0310906908696548950
    be = 0.06672455060314922
    third = 1.0/3.0
    pi34 = 0.6203504908994
    xkf = 1.919158292677513
    xks = 1.128379167095513
    # pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)

    rs = pi34/rho^third
    ec, vc = cu_XC_c_pw( rho )
    
    kf = xkf/rs
    ks = xks * CUDAnative.sqrt(kf)
    t = CUDAnative.sqrt(grho) / (2.0 * ks * rho)
    
    expe = CUDAnative.exp(-ec/ga)
    af = be / ga * (1.0 / (expe - 1.0) )
    
    bf = expe * (vc - ec)
  
    y = af * t * t
  
    xy = (1.0 + y) / (1.0 + y + y * y)
  
    qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
  
    s1 = 1.0 + be / ga * t * t * xy
    h0 = ga * CUDAnative.log(s1)
  
    dh0 = be * t * t / s1 * ( -7.0 / 3.0 * xy - qy * (af * bf / be - 7.0 / 3.0) )
  
    ddh0 = be/(2.0 * ks * ks * rho) * (xy - qy) / s1
  
    sc = rho * h0
    v1c = h0 + dh0
    v2c = ddh0

    return sc, v1c, v2c

end

function cu_XC_c_pbe_spin( rho, zeta, grho )
    # PBE correlation (without LDA part) - spin-polarized
    # iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  
    ga = 0.031091
    be = 0.06672455060314922
    third = 1.0/3.0
    pi34 = 0.6203504908994
  
    xkf = 1.919158292677513
    xks = 1.128379167095513

    # pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)


    rs = pi34/rho^third
    
    ec, vcup, vcdw = cu_XC_c_pw_spin( rho, zeta )

    kf = xkf / rs
    ks = xks * CUDAnative.sqrt(kf)
    fz = 0.5 * ( (1.0 + zeta)^(2.0 / 3.0) + (1.0 - zeta)^(2.0 / 3.0) )
    
    fz2 = fz * fz
    fz3 = fz2 * fz
    fz4 = fz3 * fz
    dfz = ( (1.0 + zeta)^(-1.0/3.0) - (1.0 - zeta)^( -1.0 / 3.0) ) / 3.0
    
    t = CUDAnative.sqrt(grho) / (2.0 * fz * ks * rho)
    expe = CUDAnative.exp( - ec / (fz3 * ga) )
    af = be / ga * (1.0 / (expe-1.0) )
    
    bfup = expe * (vcup - ec) / fz3
    bfdw = expe * (vcdw - ec) / fz3
    
    y = af * t * t
    
    xy = (1.0 + y) / (1.0 + y + y * y)
    qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
    
    s1 = 1.0 + be / ga * t * t * xy
    h0 = fz3 * ga * CUDAnative.log(s1)
    
    dh0up = be * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * (af * bfup / be - 7.0 / 3.0) )
    
    dh0dw = be * t * t * fz3 / s1 * ( - 7.0 / 3.0 * xy - qy * (af * bfdw / be - 7.0 / 3.0) )
    
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

function cu_XC_c_pw( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    a = 0.031091
    a1 = 0.21370
    b1 = 7.5957
    b2 = 3.5876
    b3 = 1.6382
    b4 = 0.49294

    # interpolation formula
    rs12 = CUDAnative.sqrt(rs)
    rs32 = rs * rs12
    rs2 = rs^2
    
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
    dom = 2.0 * a * (0.5 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 + 2.0 * b4 * rs2)
    olog = CUDAnative.log(1.0 + 1.0/om)
    ec = -2.0 * a * (1.0 + a1 * rs) * olog
    vc = -2.0*a*(1.0 + 2.0/3.0 * a1 * rs) * olog - 2.0/3.0 * a * (1.0 + a1*rs) * dom/ (om * (om + 1.0) )


    return ec, vc
end

function cu_XC_c_pw_spin( Rhoe, zeta )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    # J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
    # xc parameters, unpolarised
    a = 0.031091
    a1 = 0.21370
    b1 = 7.5957
    b2 = 3.5876
    b3 = 1.6382
    b4 = 0.49294
    c0 = a
    c1 = 0.046644
    c2 = 0.00664
    c3 = 0.01043
    d0 = 0.4335
    d1 = 1.4408
    
    # xc parameters, polarised
    ap = 0.015545
    a1p = 0.20548
    b1p = 14.1189
    b2p = 6.1977
    b3p = 3.3662
    b4p = 0.62517
    c0p = ap
    c1p = 0.025599
    c2p = 0.00319
    c3p = 0.00384
    d0p = 0.3287
    d1p = 1.7697
  
    # xc parameters, antiferro
    aa = 0.016887
    a1a = 0.11125
    b1a = 10.357
    b2a = 3.6231
    b3a = 0.88026
    b4a = 0.49671
    c0a = aa
    c1a = 0.035475
    c2a = 0.00188
    c3a = 0.00521
    d0a = 0.2240
    d1a = 0.3969
  
    fz0 = 1.709921

    #
    #     if(rs.lt.0.5d0) then
    # high density formula (not implemented)
    #
    #     else if(rs.gt.100.d0) then
    # low density formula  (not implemented)
    #
    #     else
    # interpolation formula
  
    zeta2 = zeta * zeta
    zeta3 = zeta2 * zeta
    zeta4 = zeta3 * zeta
  
    rs12 = CUDAnative.sqrt(rs)
    rs32 = rs * rs12
    rs2 = rs^2
  
    # unpolarised
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
    dom = 2.0 * a * (0.5 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 + 2.0 * b4 * rs2)
    olog = CUDAnative.log(1.0 + 1.0 / om)
    epwc = -2.0 * a * (1.0 + a1 * rs) * olog
    vpwc = -2.0 * a * (1.0 + 2.0 / 3.0 * a1 * rs) * olog - 2.0/3.0 * a * (1.0 + a1 * rs) * dom / (om * (om + 1.0) )
  
    # polarized
    omp = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
    domp = 2.0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5 * b3p * rs32 + 2.0 * b4p * rs2)
    ologp = CUDAnative.log(1.0 + 1.0 / omp)
    epwcp = -2.0 * ap * (1.0 + a1p * rs) * ologp
    vpwcp = -2.0 * ap * (1.0 + 2.0 / 3.0 * a1p * rs) * ologp - 2.0/3.0 * ap * (1.0 + a1p * rs) * domp / (omp * (omp + 1.0) )
  
    # antiferro
    oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
    doma = 2.0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5 * b3a * rs32 + 2.0 * b4a * rs2)
    ologa = CUDAnative.log(1.0 + 1.0 / oma)
    alpha = 2.0 * aa * (1.0 + a1a * rs) * ologa
    vpwca = 2.0 * aa * (1.0 + 2.0 / 3.0 * a1a * rs) * ologa + 2.0/3.0 * aa * (1.0 + a1a * rs) * doma / (oma * (oma + 1.0) )
  
    fz = ( (1.0 + zeta)^(4.0/3.0) + (1.0 - zeta)^(4.0/3.0) - 2.0) / (2.0^(4.0/3.0) - 2.0)
    dfz = ( (1.0 + zeta)^(1.0/3.0) - (1.0 - zeta)^(1.0/3.0) ) * 4.0 / (3.0 * (2.0^(4.0/3.0) - 2.0) )
  
    ec = epwc + alpha * fz * (1.0 - zeta4) / fz0 + (epwcp - epwc) * fz * zeta4
  
    vcup = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 + (vpwcp - vpwc)*fz*zeta4 +
           (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0*fz*zeta3) +
           (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3) ) * (1.0 - zeta)

    vcdw = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 + (vpwcp - vpwc) * fz * zeta4 -
          (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0 * fz * zeta3) +
            (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3) ) * (1.0 + zeta)

    return ec, vcup, vcdw
end


function cu_XC_c_vwn( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    a = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    q = CUDAnative.sqrt(4.0 * c - b * b)
    f1 = 2.0 * b / q
    f2 = b * x0 / (x0 * x0 + b * x0 + c)
    f3 = 2.0 * (2.0 * x0 + b) / q
    rs12 = CUDAnative.sqrt(rs)
    fx = rs + b * rs12 + c
    qx = CUDAnative.atan( q/(2.0 * rs12 + b) )
    ec = a * ( CUDAnative.log(rs/fx) + f1 * qx - f2 * ( CUDAnative.log( (rs12 - x0)^2 /fx) + f3 * qx) )
    
    tx = 2.0 * rs12 + b
    tt = tx * tx + q * q
    vc = ec - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b / tt -
         f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b) / tt) )

    return ec, vc

end


function cu_XC_c_vwn_spin( Rhoe, zeta )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    A      = ( 0.0310907, 0.01554535, -0.01688686394039 )
    x0     = ( -0.10498, -0.32500, -0.0047584 )
    b      = ( 3.72744, 7.06042, 1.13107 )
    c      = ( 12.9352, 18.0578, 13.0045 )
    Q      = ( 6.15199081975908, 4.73092690956011, 7.12310891781812 )
    tbQ    = ( 1.21178334272806, 2.98479352354082, 0.31757762321188 )
    fx0    = ( 12.5549141492, 15.8687885, 12.99914055888256 )
    bx0fx0 = ( -0.03116760867894, -0.14460061018521, -0.00041403379428 )

    # N.B.: A is expressed in Hartree
    # Q = sqrt(4*c - b^2)
    # tbQ = 2*b/Q
    # fx0 = X(x_0) = x_0^2 + b*x_0 + c
    # bx0fx0 = b*x_0/X(x_0)

    cfz = 2.0^(4.0/3.0) - 2.0
    cfz1 = 1.0 / cfz
    cfz2 = 4.0/3.0 * cfz1
    iddfz0 = 9.0 / 8.0 *cfz
    sqrtrs = CUDAnative.sqrt(rs)
    zeta3 = zeta^3
    zeta4 = zeta3*zeta
    trup = 1.0 + zeta
    trdw = 1.0 - zeta
    trup13 = trup^(1.0/3.0)
    trdw13 = trdw^(1.0/3.0)
    fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0)         # f(zeta)
    dfz = cfz2 * (trup13 - trdw13)     # d f / d zeta

    ecP, vcP = cu_padefit( sqrtrs, 1, x0, Q, b, c, A, tbQ, bx0fx0 )    # ecF = e_c Paramagnetic
    ecF, vcF = cu_padefit( sqrtrs, 2, x0, Q, b, c, A, tbQ, bx0fx0 )    # ecP = e_c Ferromagnetic
    ac, dac  = cu_padefit( sqrtrs, 3, x0, Q, b, c, A, tbQ, bx0fx0 )    # ac = "spin stiffness"

    ac = ac * iddfz0
    dac = dac * iddfz0
    De = ecF - ecP - ac # e_c[F] - e_c[P] - alpha_c/(ddf/ddz(z=0))
    fzz4 = fz * zeta4
    ec = ecP + ac * fz  + De * fzz4

    dec1 = vcP + dac*fz + (vcF - vcP - dac) * fzz4       # e_c - (r_s/3)*(de_c/dr_s)
    dec2 = ac*dfz + De*(4.0*zeta3*fz + zeta4*dfz)      # de_c/dzeta

    # v_c[s] = e_c - (r_s/3)*(de_c/dr_s) + [sign(s)-zeta]*(de_c/dzeta)
    vcup = dec1 + (1.0 - zeta)*dec2
    vcdw = dec1 - (1.0 + zeta)*dec2

    return ec, vcup, vcdw
end


function cu_padefit(x, i, x0, Q, b, c, A, tbQ, bx0fx0)
    # implements formula [4.4] in:
    # S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)

    # Pade fit calculated in x and its derivative w.r.t. rho
    # rs = inv((rho*)^(1/3)) = x^2
    # fit  [eq. 4.4]
    # dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx

    sqx = x * x                          # x^2 = r_s
    xx0 = x - x0[i]                      # x - x_0
    Qtxb = Q[i] / (2.0*x + b[i])      # Q / (2x+b)
    atg = CUDAnative.atan(Qtxb)                     # tan^-1(Q/(2x+b))
    fx = sqx + b[i]*x + c[i]             # X(x) = x^2 + b*x + c

    fit = A[i] * (  CUDAnative.log(sqx/fx) + tbQ[i]*atg - 
            bx0fx0[i] * ( CUDAnative.log(xx0*xx0/fx) + (tbQ[i] + 4.0*x0[i]/Q[i]) * atg )  )

    txb = 2.0*x + b[i]
    txbfx = txb / fx
    itxbQ = 1.0 / (txb*txb + Q[i]*Q[i])

    dfit = fit - A[i] / 3.0 + A[i]*x/6.0 * (  txbfx + 4.0*b[i]*itxbQ + 
              bx0fx0[i] * ( 2.0/xx0 - txbfx - 4.0*( b[i] + 2.0*x0[i])*itxbQ )  )


    return fit, dfit
end


# PBE exchange (without Slater exchange):
# iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
function cu_XC_x_pbe( rho, grho )

    third = 1.0/3.0
    k = 0.804
    mu = 0.2195149727645171
    c1 = 0.75/pi
    c2 = 3.093667726280136
    c5 = 4.0*third

    agrho = CUDAnative.sqrt(grho)
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
    sx = exunif * fx

    # Potential
    dxunif = exunif * third
    dfx1 = f2 * f2
    dfx = 2.0 * mu * s1 / dfx1
    
    v1x = sx + dxunif * fx + exunif * dfx * ds
    v2x = exunif * dfx * dsg / agrho
    sx = sx * rho

    return sx, v1x, v2x
end


function cu_XC_x_slater( Rhoe::Float64 )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    f = -0.687247939924714
    alpha = 2.0/3.0

    ex = f * alpha / rs
    vx = 4.0 / 3.0 * f * alpha / rs
    return ex, vx
end


function cu_XC_x_slater_spin( rho, zeta )
    #     Slater exchange with alpha=2/3, spin-polarized case
    #
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
