# PBE correlation (without LDA part)
# iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).

# rewritten to follow closely the paper

function XC_c_pbe_v2( ρ, ∇ρ2 )

    γ = (1 - log(2))/pi^2
    β = 0.06672455060314922
    β_div_γ = β/γ
    THIRD = 1.0/3.0

    rs = ( 3.0/(4.0*pi*ρ) )^THIRD
    ec, vc = XC_c_pw( ρ )
    
    kf = (9.0*pi/4.0)^THIRD / rs
    ks = sqrt(4.0*kf/pi)

    t = sqrt(∇ρ2) / (2.0*ks*ρ)
    A = β_div_γ * (1.0 / (exp(-ec/γ) - 1.0) )
    At2 = A*t^2
    H = γ * log(1 + β_div_γ * t^2 * (1 + At2)/(1 + At2 + At2^2) )
  
    #bf = exp(-ec/γ) * (vc - ec)
    #y = A * t * t
    #xy = (1.0 + y) / (1.0 + y + y * y)
    #qy = y * y * (2.0 + y) / (1.0 + y + y * y)^2
    #s1 = 1.0 + β / γ * t^2 * xy

    #dH = β * t * t / s1 * ( -7.0 / 3.0 * xy - qy * (A * bf / β - 7.0 / 3.0) )
  
    #ddH = β/(2.0 * ks * ks * ρ) * (xy - qy) / s1
  
    sc = ρ * H
    
    v1c = 0.0 #H + dH
    
    v2c = 0.0 #ddH

    return sc, v1c, v2c

end

