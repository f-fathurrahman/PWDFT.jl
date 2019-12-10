function kernel_eval_Vloc_R_GTH( Zval, rlocal, C1, C2, C3, C4,
    r::CuArray{Float64,1},
    Vloc::CuArray{Float64,1}
)
    ip = ( blockIdx().x - 1 ) * blockDim().x + threadIdx().x
    if ip <= length(r)
        rrloc = r[ip]/rlocal
        term1 = C1 + C2*rrloc^2 + C3*rrloc^4 + C4*rrloc^6
        Vloc[ip] = -Zval/r[ip] * CUDAnative.erf( rrloc/sqrt(2.0) ) + exp(-0.5*rrloc^2)*term1
    end
    return
end


function kernel_eval_Vloc_G_GTH( Zval, rlocal, C1, C2, C3, C4,
    G2::CuArray{Float64,1},
    Vg::CuArray{Float64,1}
)

    pre1 = -4*pi*Zval
    pre2 = CUDAnative.sqrt(8*pi^3)*rlocal^3

    ig = ( blockIdx().x - 1 ) * blockDim().x + threadIdx().x

    #
    # VG[1] = V(G2=0) should be set to zero
    #
    if (ig != 1) && ig >= length(G2)
        Gr = CUDAnative.sqrt( G2[ig] )*rlocal
        expGr2 = CUDAnative.exp( -0.5*Gr^2 )
        Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (C1 + C2*(3.0 - Gr^2) +
                 C3*(15.0 - 10.0*Gr^2 + Gr^4) + C4*(105.0 - 105.0*Gr^2 + 21.0*Gr^4 - Gr^6) )
    end

    return
end

# rrl = psp.rc[l+1]
function kernel_eval_proj_G_GTH(
    rrl::Float64, l::Int64, iproj::Int64, CellVolume::Float64,
    Gm::CuArray{Float64,1},
    Vprj::CuArray{Float64,1}
)

    ig = ( blockIdx().x - 1 ) * blockDim().x + threadIdx().x

    perOmega = CUDAnative.pow(2.0, (l+1)) * CUDAnative.pow(rrl, (2*l+3)) / CellVolume
    pre =  4.0 * CUDAnative.pow(Float64(pi), (5.0/4.0)) * CUDAnative.sqrt( perOmega )

    Ng = length(Gm)

    # s-channel
    if l == 0

        if iproj==1
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * CUDAnative.exp( -0.5*Gr2 )
            end

        elseif iproj==2
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * 2.0/CUDAnative.sqrt(15.0) * CUDAnative.exp( -0.5*Gr2 ) * ( 3.0 - Gr2 )
            end

        elseif iproj==3
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (4.0/3.0)/CUDAnative.sqrt(105.0) * CUDAnative.exp( -0.5*Gr2 ) * (15.0 - 10.0*Gr2 + Gr2^2)
            end

        end  # if iproj

    # p-channel
    elseif l == 1

        if iproj == 1
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (1.0/CUDAnative.sqrt(3.0)) * CUDAnative.exp(-0.5*Gr2) * Gm[ig]
            end

        elseif iproj == 2
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (2.0/CUDAnative.sqrt(105.0)) * CUDAnative.exp(-0.5*Gr2) * Gm[ig] * (5.0 - Gr2)
            end

        elseif iproj == 3
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (4.0/3.0)/CUDAnative.sqrt(1155.0) * CUDAnative.exp(-0.5*Gr2) * Gm[ig] * (35.0 - 14.0*Gr2 + Gr2^2)
            end

        end # if iproj

    # d-channel
    elseif l == 2

        if iproj == 1
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (1.0/CUDAnative.sqrt(15.0)) * CUDAnative.exp(-0.5*Gr2) * Gm[ig]^2
            end

        elseif iproj == 2
            if ig <= Ng
                Gr2 = ( Gm[ig]*rrl )^2
                Vprj[ig] = pre * (2.0/3.0)/CUDAnative.sqrt(105.0) * CUDAnative.exp(-0.5*Gr2) * Gm[ig]^2 * (7.0 - Gr2)
            end

        end # if iproj

    # f-channel
    elseif l == 3
        # XXX only one projector
        if ig <= Ng
            Gr2 = ( Gm[ig]*rrl )^2
            Vprj[ig] = pre * Gm[ig]^3 * CUDAnative.exp(-0.5*Gr2)/CUDAnative.sqrt(105.0)
        end

    end  # if l

    return
end


# This is meant to called from device (GPU)
function cu_eval_proj_G( rrl, l::Int64, iproj::Int64, Gm::Float64, CellVolume::Float64 )

    Vprj = 0.0

    # s-channel
    if l == 0
        if iproj==1
            Gr2 = ( Gm*rrl )^2
            Vprj = CUDAnative.exp( -0.5*Gr2 )
        elseif iproj==2
            Gr2 = ( Gm*rrl )^2
            Vprj = 2.0/CUDAnative.sqrt(15.0) * CUDAnative.exp( -0.5*Gr2 ) * ( 3.0 - Gr2 )
        elseif iproj==3
            Gr2 = ( Gm*rrl )^2
            Vprj = (4.0/3.0)/CUDAnative.sqrt(105.0) * CUDAnative.exp( -0.5*Gr2 ) * (15.0 - 10.0*Gr2 + Gr2^2)
        end  # if iproj

    # p-channel
    elseif l == 1
        if iproj == 1
            Gr2 = ( Gm*rrl )^2
            Vprj = (1.0/CUDAnative.sqrt(3.0)) * CUDAnative.exp(-0.5*Gr2) * Gm
        elseif iproj == 2
            Gr2 = ( Gm*rrl )^2
            Vprj = (2.0/CUDAnative.sqrt(105.0)) * CUDAnative.exp(-0.5*Gr2) * Gm * (5.0 - Gr2)
        elseif iproj == 3
            Gr2 = (  Gm*rrl )^2
            Vprj = (4.0/3.0)/CUDAnative.sqrt(1155.0) * CUDAnative.exp(-0.5*Gr2) * Gm * (35.0 - 14.0*Gr2 + Gr2^2)
        end # if iproj

    # d-channel
    #elseif l == 2
    #    if iproj == 1
    #        Gr2 = ( Gm*rrl )^2
    #        Vprj = (1.0/CUDAnative.sqrt(15.0)) * CUDAnative.exp(-0.5*Gr2) * Gm^2
    #    elseif iproj == 2
    #        Gr2 = ( Gm*rrl )^2
    #        Vprj = (2.0/3.0)/CUDAnative.sqrt(105.0) * CUDAnative.exp(-0.5*Gr2) * Gm^2 * (7.0 - Gr2)
    #    end # if iproj

    # f-channel
    #elseif l == 3
    #    # XXX only one projector
    #    Gr2 = ( Gm*rrl )^2
    #    Vprj = Gm^3 * CUDAnative.exp(-0.5*Gr2)/CUDAnative.sqrt(105.0)

    end  # if l

    perOmega = CUDAnative.pow(2.0, (l+1)) * CUDAnative.pow(rrl, (2*l+3)) / CellVolume
    pre =  4.0 * CUDAnative.pow(Float64(pi), (5.0/4.0)) * CUDAnative.sqrt( perOmega )
    Vprj = pre * Vprj
    return Vprj
end
