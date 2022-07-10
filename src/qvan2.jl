# TODO: using pspotNL
function qvan2!(
    pspotNL,
    ih, jh, isp, G2, ylmk0, QfuncG
)

    qvan2!(
        pspotNL.indv, pspotNL.nhtolm, pspotNL.lpl, pspotNL.lpx, pspotNL.ap,
        pspotNL.qradG,
        ih, jh, isp, G2, ylmk0, QfuncG
    )
    
    return
end




function qvan2!(
    indv, nhtolm, lpl, lpx, ap,
    qradG,
    ih, jh, isp, G2, ylmk0, QfuncG
)
    Ng = size(G2,1)

    nb = indv[ih,isp]
    mb = indv[jh,isp]

    if nb >= mb
        ijv = round(Int64, nb * (nb - 1) / 2 + mb)
    else
        ijv = round(Int64, mb * (mb - 1) / 2 + nb)
    end

    ivl = nhtolm[ih,isp]
    jvl = nhtolm[jh,isp]

    dq = 0.01 # XXX HARDCODED

    fill!(QfuncG, 0.0) # zero out QfuncG

    for lm in 1:lpx[ivl,jvl]
        #
        lp = lpl[ivl,jvl,lm] # combined (l+1) index (using default 1-index based array)

        # Hardcoded assertion
        @assert lp >= 1
        @assert lp <= 49

        # finds angular momentum l corresponding to combined index lp (l is 
        # actually l+1 because this is the way qrad is stored, check init_us_1)
        if lp == 1
            l = 1
        elseif lp <= 4
            l = 2
        elseif lp <= 9
            l = 3
        elseif lp <= 16
            l = 4
        elseif lp <= 25
            l = 5
        elseif lp <= 36
            l = 6
        else
            l = 7
        end
        
        # sig = (-im)^(l-1)  # using physics' l (start from 0)
        prefact = (-im)^(l-1) * ap[lp,ivl,jvl]

        for ig in 1:Ng
            # calculates quantites depending on the module of G only when needed
            Gm = sqrt(G2[ig])
            px = Gm/dq - floor(Int64, Gm/dq)
            ux = 1.0 - px
            vx = 2.0 - px
            wx = 3.0 - px
            i0 = floor(Int64, Gm/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
            uvx = ux * vx * (1/6)
            pwx = px * wx * 0.5
            work = qradG[isp][i0,ijv,l] * uvx * wx +
                   qradG[isp][i1,ijv,l] * pwx * vx -
                   qradG[isp][i2,ijv,l] * pwx * ux +
                   qradG[isp][i3,ijv,l] * px * uvx
            QfuncG[ig] = QfuncG[ig] + prefact * ylmk0[ig,lp] * work
        end
    end

    return
end