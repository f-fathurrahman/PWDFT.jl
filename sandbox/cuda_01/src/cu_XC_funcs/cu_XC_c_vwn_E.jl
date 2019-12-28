function cu_XC_c_vwn_E( Rhoe )

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)
    rs = pi34/CUDAnative.pow(Rhoe, third)

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

    return ec

end