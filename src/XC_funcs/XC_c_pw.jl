function XC_c_pw( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994
    rs = pi34/Rhoe^third

    a = 0.031091
    a1 = 0.21370
    b1 = 7.5957
    b2 = 3.5876
    b3 = 1.6382
    b4 = 0.49294

    # interpolation formula
    rs12 = sqrt(rs)
    rs32 = rs * rs12
    rs2 = rs^2
    
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
    dom = 2.0 * a * (0.5 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 + 2.0 * b4 * rs2)
    olog = log(1.0 + 1.0/om)
    ec = -2.0 * a * (1.0 + a1 * rs) * olog
    vc = -2.0*a*(1.0 + 2.0/3.0 * a1 * rs) * olog - 2.0/3.0 * a * (1.0 + a1*rs) * dom/ (om * (om + 1.0) )

    return ec, vc
end