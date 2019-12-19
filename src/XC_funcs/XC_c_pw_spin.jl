function XC_c_pw_spin( Rhoe, zeta )

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
  
    rs12 = sqrt(rs)
    rs32 = rs * rs12
    rs2 = rs^2
  
    # unpolarised
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
    dom = 2.0 * a * (0.5 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 + 2.0 * b4 * rs2)
    olog = log(1.0 + 1.0 / om)
    epwc = -2.0 * a * (1.0 + a1 * rs) * olog
    vpwc = -2.0 * a * (1.0 + 2.0 / 3.0 * a1 * rs) * olog - 2.0/3.0 * a * (1.0 + a1 * rs) * dom / (om * (om + 1.0) )
  
    # polarized
    omp = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
    domp = 2.0 * ap * (0.5 * b1p * rs12 + b2p * rs + 1.5 * b3p * rs32 + 2.0 * b4p * rs2)
    ologp = log(1.0 + 1.0 / omp)
    epwcp = -2.0 * ap * (1.0 + a1p * rs) * ologp
    vpwcp = -2.0 * ap * (1.0 + 2.0 / 3.0 * a1p * rs) * ologp - 2.0/3.0 * ap * (1.0 + a1p * rs) * domp / (omp * (omp + 1.0) )
  
    # antiferro
    oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
    doma = 2.0 * aa * (0.5 * b1a * rs12 + b2a * rs + 1.5 * b3a * rs32 + 2.0 * b4a * rs2)
    ologa = log(1.0 + 1.0 / oma)
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