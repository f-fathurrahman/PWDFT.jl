function calc_exx_divergence(pw,
    Nq1, Nq2, Nq3,
    x_gamma_extrapolation,
    erf_scrlen,
    erfc_scrlen,
    yukawa,
    grid_factor,
    nqs;
    use_regularization = true
)
    
    if !use_regularization
        return 0.0
    end

    ecutwfc = pw.ecutwfc
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    Ng = pw.gvec.Ng
    G = pw.gvec.G
    CellVolume = pw.CellVolume

    alat = sqrt(LatVecs[1,1]^2 + LatVecs[2,1]^2 + LatVecs[3,1]^2)
    tpiba2  = (2π/alat)^2

    gcutw = 2*ecutwfc/tpiba2
    α = 10.0/gcutw

    #println("alat = ", alat)
    #println("tpiba2 = ", tpiba2)
    #println("gcutw = ", gcutw)
    #println("α = ", α)

    SMALL = 1e-6

    dq1 = 1.0/Nq1
    dq2 = 1.0/Nq2 
    dq3 = 1.0/Nq3 
    res = 0.0
    
    odg = zeros(Bool, 3)
    xq = zeros(Float64, 3)
    q = zeros(Float64, 3)

    on_double_grid = false
    nqhalf = 0.5*[Nq1, Nq2, Nq3]

    for iq1 in 1:Nq1, iq2 in 1:Nq2, iq3 in 1:Nq3
        xq[:] = RecVecs[:,1] * (iq1-1) * dq1 +
                RecVecs[:,2] * (iq2-1) * dq2 +
                RecVecs[:,3] * (iq3-1) * dq3
        for ig in 1:Ng
            q[1] = xq[1] + G[1,ig]
            q[2] = xq[2] + G[2,ig]
            q[3] = xq[3] + G[3,ig]
            qq = (q[1]^2 + q[2]^2 + q[3]^2)/tpiba2 # XXX scale this ???
            if x_gamma_extrapolation
                x = ( q[1]*LatVecs[1,1] + q[2]*LatVecs[2,1] + q[3]*LatVecs[3,1] ) * nqhalf[1] / (2*pi)
                odg[1] = abs(x - round(Int64,x)) < SMALL
                #
                x = ( q[1]*LatVecs[1,2] + q[2]*LatVecs[2,2] + q[3]*LatVecs[3,2] ) * nqhalf[2] / (2*pi)
                odg[2] = abs(x - round(Int64,x)) < SMALL
                #
                x = ( q[1]*LatVecs[1,3] + q[2]*LatVecs[2,3] + q[3]*LatVecs[3,3] ) * nqhalf[3] / (2*pi)
                odg[3] = abs(x - round(Int64,x)) < SMALL
                #
                on_double_grid = all(odg)
            end
            if !on_double_grid
                if qq > 1e-8 # XXX Note that qq is in tpiba2 unit
                    if erfc_scrlen > 0
                        res += exp(-α*qq) / qq * (1.0 - exp(-qq/4.0/erfc_scrlen^2)) * grid_factor
                    elseif erf_scrlen > 0
                        res += exp(-α*qq) / qq * (exp(-qq/4.0/erf_scrlen^2)) * grid_factor
                    else
                        #XXX qq is in tpiba2 unit, so `yukawa` also need to be converted
                        res += exp(-α*qq) / (qq + yukawa/tpiba2) * grid_factor
                    end
                    #println("$ig $qq $res")
                end
            end
        end
    end
    #println("grid_factor = ", grid_factor)
    #println("Line 86 res = ", res)
    
    if !x_gamma_extrapolation
        if yukawa > 0.0
            res += 4*pi^2/yukawa
        elseif erfc_scrlen > 0.0
            res += pi^2/erfc_scrlen^2
        else
            res -= α
        end
    end
    #println("Line 97 res = ", res)

    res = res*alat^2 / pi
    #res = res*pi/nqs
    #println("Line 111: res (in Ry) = ", 2*res)

    #α = α / (4*pi^2)
    α = α / tpiba2
    nqq = 100000
    dq = 5.0 / sqrt(α) / nqq
    aa = 0.0
    for iq in 0:nqq
        q_ = dq * (iq + 0.5)
        qq = q_ * q_ / tpiba2 # XXX
        if erfc_scrlen > 0
            aa -= exp(-α*qq) * exp(-qq/4/erfc_scrlen^2)*dq
        elseif erf_scrlen > 0
            aa = 0.0
        else
            aa -= -exp(-α*qq) * yukawa / (qq + yukawa)*dq
        end
    end
    #println("Line 129: aa = ", aa)
    aa *= 8/(4*pi)
    aa += 1.0/sqrt(α*π)
    #println("Line 132 aa = ", aa)
    if erf_scrlen > 0
        aa = 1.0/sqrt( (α + 1.0/4.0/erf_scrlen^2) * pi )
    end
    res -= CellVolume*aa
    #println("res in Ry = ", res*2)

    return res*nqs
end