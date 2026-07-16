function g2_convolution!( exx, LatVecs, xk, xkq, fac )

    Ng = exx.gvec.Ng
    G = exx.gvec.G
    SMALL = 1e-6
    SMALL_QDIV = 1e-8 # from exx
    grid_factor = exx.grid_factor

    gau_scrlen = exx.gau_scrlen
    erfc_scrlen = exx.erfc_scrlen
    erf_scrlen = exx.erf_scrlen
    x_gamma_extrapolation = exx.x_gamma_extrapolation
    exxdiv = exx.exxdiv
    yukawa = exx.yukawa
    Nq1 = exx.Nq1
    Nq2 = exx.Nq2
    Nq3 = exx.Nq3

    q = zeros(Float64, 3)
    grid_factor_track = zeros(Float64, Ng)
    qq_track = zeros(Float64, Ng)
    odg = zeros(Bool, 3)

    alat = sqrt(LatVecs[1,1]^2 + LatVecs[2,1]^2 + LatVecs[3,1]^2)
    tpiba2  = (2π/alat)^2

    # First the types of Coulomb potential that need q(3) and an external call
    #=
    if use_coulomb_vcut_ws
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            fac[ig] = vcut_get(vcut,q)
        end
        return
    end
    if use_coulomb_vcut_spheric
        for ig in 1:Ng
            q[:]= xk[:] - xkq[:] + G[:,ig]
            fac(ig) = vcut_spheric_get(vcut,q)
        end
        return
    end
    =#

    #
    # Now the Coulomb potential that are computed on the fly
    nqhalf = 0.5*[Nq1, Nq2, Nq3]
    #
    # Set the grid_factor_track and qq_track
    #
    #XXX: Instead of LatVecs we actually need inv(RecVecs)
    if x_gamma_extrapolation
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            qq_track[ig] = sum(q.^2)
            #
            x1 = ( q[1]*LatVecs[1,1] + q[2]*LatVecs[2,1] + q[3]*LatVecs[3,1] ) * nqhalf[1] / (2*pi)
            odg[1] = abs(x1 - round(Int64, x1)) < SMALL
            #
            x2 = ( q[1]*LatVecs[1,2] + q[2]*LatVecs[2,2] + q[3]*LatVecs[3,2] ) * nqhalf[2] / (2*pi)
            odg[2] = abs(x2 - round(Int64, x2)) < SMALL
            #
            x3 = ( q[1]*LatVecs[1,3] + q[2]*LatVecs[2,3] + q[3]*LatVecs[3,3] ) * nqhalf[3] / (2*pi)
            odg[3] = abs(x3 - round(Int64, x3)) < SMALL
            #
            if all(odg)
                grid_factor_track[ig] = 0.0 # on double grid
            else
                grid_factor_track[ig] = grid_factor # not on double grid
            end
            #println("x1=$x1 x2=$x2 x3=$x3 odg=$odg")
        end
    else
        # No gamma extrapolation
        for ig in 1:Ng
            q[:] = xk[:] - xkq[:] + G[:,ig]
            qq_track[ig] = sum(q.^2)
        end
        fill!(grid_factor_track, 1.0)
    end
    #println("xk = ", xk)
    #println("xkq = ", xkq)
    #println("sum qq_track = ", sum(qq_track))
    #println("sum grid_factor_track = ", sum(grid_factor_track))
    #
    # The big loop
    for ig in 1:Ng
        qq = qq_track[ig]
        #
        if gau_scrlen > 0
            fac[ig] = ( (pi/gau_scrlen)^1.5 )*exp(-qq/4/gau_scrlen) * grid_factor_track[ig]
        #
        elseif qq/tpiba2 > SMALL_QDIV #XXX eps_qdiv is given in tpiba2 ?
            if erfc_scrlen > 0
                fac[ig] = 4*pi/qq*(1.0 - exp(-qq/4/erfc_scrlen^2)) * grid_factor_track[ig]
            elseif erf_scrlen > 0
                fac[ig] = 4*pi/qq*( exp(-qq/4/erf_scrlen^2) ) * grid_factor_track[ig]
            else
                fac[ig] = 4*pi/( qq + yukawa ) * grid_factor_track[ig]
            end
        #
        else
            # Small q
            #
            fac[ig] = -exxdiv # or rather something else (see F.Gygi)
            if yukawa > 0.0 && !x_gamma_extrapolation
                fac[ig] += 4*pi/( qq + yukawa ) # add additional factor
            end
            #      
            if erfc_scrlen > 0.0 && !x_gamma_extrapolation
                fac[ig] += pi/(erfc_scrlen^2)
            end
        end
    end

    return
end