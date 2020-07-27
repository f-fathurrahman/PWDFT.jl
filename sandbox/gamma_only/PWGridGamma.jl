function inv_mm_to_nn(nn::Int64, S::Int64)
    if nn < 0
        return nn + S
    else
        return nn
    end
end


struct GVectorsGamma
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    idx_g2rm::Array{Int64,1}
    G2_shells::Array{Float64,1}
    idx_g2shells::Array{Int64,1}
end


#
# Based on ggen of QE-6.5
#

function GVectorsGamma( Ns, RecVecs, ecutrho )

    ni = floor( Int64, (Ns[1]-1)/2 )
    nj = floor( Int64, (Ns[2]-1)/2 )
    nk = floor( Int64, (Ns[3]-1)/2 )
    
    # gamma-only: exclude space with x < 0
    istart = 0

    G_tmp = zeros(Float64,3)
    tt = zeros(Float64,Ns[3])

    ig = 0
    ip = 0
    ipm = 0

    SMALL = eps()

    Ng = PWDFT.calc_Ng( Ns, RecVecs, ecutrho )
    Ng = round(Int64, (Ng + 1)/2)

    G  = zeros(Float64,3,Ng)
    G2 = zeros(Float64,Ng)
    idx_g2r = zeros(Int64,Ng)  # G=0 is included here
    idx_g2rm = zeros(Int64,Ng) # for negative of G, NOTE: ig=1 should not be accessed

    for i in istart:ni

        if i == 0
           jstart = 0
        else
           jstart = -nj
        end
       
        for j in jstart:nj

            if ( (i == 0) && (j == 0) )
                kstart = 0
            else
                kstart = -nk
            end
            
            for k in kstart:nk
                
                G_tmp[1] = RecVecs[1,1]*i + RecVecs[1,2]*j + RecVecs[1,3]*k
                G_tmp[2] = RecVecs[2,1]*i + RecVecs[2,2]*j + RecVecs[2,3]*k
                G_tmp[3] = RecVecs[3,1]*i + RecVecs[3,2]*j + RecVecs[3,3]*k
                G2_tmp = G_tmp[1]^2 + G_tmp[2]^2 + G_tmp[3]^2
                
                if 0.5*G2_tmp <= ecutrho
                    #
                    ig = ig + 1
                    #
                    G[1,ig] = G_tmp[1]
                    G[2,ig] = G_tmp[2]
                    G[3,ig] = G_tmp[3]
                    #
                    G2[ig] = G2_tmp
                    #
                    ip1 = inv_mm_to_nn(i, Ns[1])
                    ip2 = inv_mm_to_nn(j, Ns[2])
                    ip3 = inv_mm_to_nn(k, Ns[3])
                    #
                    ip = ip1 + 1 + Ns[1]*( ip2 + Ns[2]*ip3 )
                    #
                    idx_g2r[ig] = ip # index of +G
                    #
                    # Index of -G
                    if G2_tmp > SMALL
                        ip1m = inv_mm_to_nn(-i, Ns[1])
                        ip2m = inv_mm_to_nn(-j, Ns[2])
                        ip3m = inv_mm_to_nn(-k, Ns[3])
                        ipm = ip1m + 1 + Ns[1]*( ip2m + Ns[2]*ip3m )
                        idx_g2rm[ig] = ipm
                    end
                end # if
            end # kstart:nk

       end
    end

    idx_sorted = sortperm(G2)
    G = G[:,idx_sorted]
    G2 = G2[idx_sorted]
    idx_g2r = idx_g2r[idx_sorted]
    idx_g2rm = idx_g2rm[idx_sorted]

    G2_shells, idx_g2shells = PWDFT.init_Gshells( G2 )

    return GVectorsGamma( Ng, G, G2, idx_g2r, idx_g2rm, G2_shells, idx_g2shells )
end

function calc_Ngw_gamma( ecutwfc::Float64, gvec::GVectorsGamma )
    G = gvec.G
    Ng = gvec.Ng
    # k = [0,0,0]
    Ngw = 0
    for ig = 1:Ng
        Gk2 = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
        if 0.5*Gk2 <= ecutwfc
            Ngw = Ngw + 1
        end
    end
    return Ngw
end

struct GVectorsWGamma
    Ngw::Int64
    idx_gw2g::Array{Int64,1}
    idx_gw2r::Array{Int64,1}  # for half G vectors set
    idx_gw2rm::Array{Int64,1} # for other half (minus) G vectors set
end

function GVectorsWGamma( ecutwfc::Float64, gvec::GVectorsGamma )
    G = gvec.G
    Ngw = calc_Ngw_gamma(ecutwfc, gvec)
    idx_gw2g = zeros(Int64,Ngw)
    idx_gw2r = zeros(Int64,Ngw)
    idx_gw2rm = zeros(Int64,Ngw)
    #
    igw = 0
    for ig = 1:Ngw
        Gk2 = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
        if 0.5*Gk2 <= ecutwfc
            igw = igw + 1
            idx_gw2g[igw] = ig
            idx_gw2r[igw] = gvec.idx_g2r[ig]
            idx_gw2rm[igw] = gvec.idx_g2rm[ig]
        end
    end
    return GVectorsWGamma( Ngw, idx_gw2g, idx_gw2r, idx_gw2rm )
end

struct PWGridGamma
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    gvec::GVectorsGamma
    gvecw::GVectorsWGamma
    planfw::PWDFT.PLANFW_TYPE
    planbw::PWDFT.PLANBW_TYPE
end


function PWGridGamma( ecutwfc::Float64, LatVecs::Array{Float64,2}; Ns_=(0,0,0) )

    ecutrho = 4.0*ecutwfc
    #
    RecVecs = 2*pi*inv(Matrix(LatVecs'))

    CellVolume = abs(det(LatVecs))
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    if any(Ns_ .== 0)
        Ns1 = good_fft_order(Ns1)
        Ns2 = good_fft_order(Ns2)
        Ns3 = good_fft_order(Ns3)
        Ns = (Ns1,Ns2,Ns3)
    else
        Ns = Ns_[:]
    end

    Npoints = prod(Ns)
    
    gvec = GVectorsGamma( Ns, RecVecs, ecutrho )

    gvecw = GVectorsWGamma( ecutwfc, gvec )

    planfw = plan_fft( zeros(ComplexF64,Ns) )
    planbw = plan_ifft( zeros(ComplexF64,Ns) )

    return PWGridGamma(
        ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, gvec, gvecw,
        planfw, planbw
    )
end