function XC_c_vwn( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/Rhoe^third

    a = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    q = sqrt(4.0 * c - b * b)
    f1 = 2.0 * b / q
    f2 = b * x0 / (x0 * x0 + b * x0 + c)
    f3 = 2.0 * (2.0 * x0 + b) / q
    rs12 = sqrt(rs)
    fx = rs + b * rs12 + c
    qx = atan( q/(2.0 * rs12 + b) )
    ec = a * ( log(rs/fx) + f1 * qx - f2 * ( log( (rs12 - x0)^2 /fx) + f3 * qx) )
    
    tx = 2.0 * rs12 + b
    tt = tx * tx + q * q
    vc = ec - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b / tt -
         f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b) / tt) )

    return ec, vc

end