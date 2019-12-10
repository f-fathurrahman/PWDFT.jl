function eval_Vloc_R( Zval, C1, C2, C3, C4, rlocal, r::Float64 )
    term1 = C1
    rrloc = r/rlocal
    term1 = C1 + C2*rrloc^2 + C3*rrloc^4 + C4*rrloc^6
    Vloc = -Zval/r * erf( rrloc/sqrt(2.0) ) + exp(-0.5*rrloc^2)*term1
    return Vloc
end
