function XC_c_vwn_spin( Rhoe, zeta )

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
    sqrtrs = sqrt(rs)
    zeta3 = zeta^3
    zeta4 = zeta3*zeta
    trup = 1.0 + zeta
    trdw = 1.0 - zeta
    trup13 = trup^(1.0/3.0)
    trdw13 = trdw^(1.0/3.0)
    fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0)         # f(zeta)
    dfz = cfz2 * (trup13 - trdw13)     # d f / d zeta

    ecP, vcP = padefit( sqrtrs, 1, x0, Q, b, c, A, tbQ, bx0fx0 )    # ecF = e_c Paramagnetic
    ecF, vcF = padefit( sqrtrs, 2, x0, Q, b, c, A, tbQ, bx0fx0 )    # ecP = e_c Ferromagnetic
    ac, dac  = padefit( sqrtrs, 3, x0, Q, b, c, A, tbQ, bx0fx0 )    # ac = "spin stiffness"

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


function padefit(x, i, x0, Q, b, c, A, tbQ, bx0fx0)
    # implements formula [4.4] in:
    # S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)

    # Pade fit calculated in x and its derivative w.r.t. rho
    # rs = inv((rho*)^(1/3)) = x^2
    # fit  [eq. 4.4]
    # dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx

    sqx = x * x                          # x^2 = r_s
    xx0 = x - x0[i]                      # x - x_0
    Qtxb = Q[i] / (2.0*x + b[i])      # Q / (2x+b)
    atg = atan(Qtxb)                     # tan^-1(Q/(2x+b))
    fx = sqx + b[i]*x + c[i]             # X(x) = x^2 + b*x + c

    fit = A[i] * (  log(sqx/fx) + tbQ[i]*atg - 
            bx0fx0[i] * ( log(xx0*xx0/fx) + (tbQ[i] + 4.0*x0[i]/Q[i]) * atg )  )

    txb = 2.0*x + b[i]
    txbfx = txb / fx
    itxbQ = 1.0 / (txb*txb + Q[i]*Q[i])

    dfit = fit - A[i] / 3.0 + A[i]*x/6.0 * (  txbfx + 4.0*b[i]*itxbQ + 
              bx0fx0[i] * ( 2.0/xx0 - txbfx - 4.0*( b[i] + 2.0*x0[i])*itxbQ )  )


    return fit, dfit
end