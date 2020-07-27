function inv_mm_to_nn(nn::Int64, S::Int64)
    if nn < 0
        return nn + S
    else
        return nn
    end
end


#
# Based on ggen of QE-6.5
#

function init_gvec_gamma( Ns, RecVecs, ecutrho )

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
    idx_g2r = zeros(Int64,Ng) # G=0 is included here
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
                    @printf("idx_miller = [%4d %4d %4d] ", i, j, k)
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
                    @printf(" idx grid = [%4d %4d %4d] ", ip1, ip2, ip3)
                    #
                    ip = ip1 + 1 + Ns[1]*( ip2 + Ns[2]*ip3 )
                    #
                    idx_g2r[ig] = ip # index of +G
                    #
                    # Index of -G
                    if G2_tmp < SMALL # printing stuffs
                        println()
                    end
                    #
                    if G2_tmp > SMALL
                        ip1m = inv_mm_to_nn(-i, Ns[1])
                        ip2m = inv_mm_to_nn(-j, Ns[2])
                        ip3m = inv_mm_to_nn(-k, Ns[3])
                        ipm = ip1m + 1 + Ns[1]*( ip2m + Ns[2]*ip3m )
                        idx_g2rm[ig] = ipm
                        @printf("idx negative = [%4d %4d %4d] ipm = %4d\n", ip1m, ip2m, ip3m, ipm)
                    end
                end # if
            end # kstart:nk

       end
    end

    println("Last ig = ", ig)
    println("Ng = ", Ng)

    println("Before sort")
    for ig in 2:6
        @printf("ig = %4d G = [%10.5f,%10.5f,%10.5f] G2 = %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    for ig in Ng-5:Ng
        @printf("ig = %4d G = [%10.5f,%10.5f,%10.5f] G2 = %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end

    idx_sorted = sortperm(G2)
    G = G[:,idx_sorted]
    G2 = G2[idx_sorted]
    idx_g2r = idx_g2r[idx_sorted]
    idx_g2rm = idx_g2rm[idx_sorted]

    println("max idx_g2r : ", maximum(idx_g2r))
    println("max idx_g2rm: ", maximum(idx_g2rm))
    println("Npoints = ", prod(Ns))

    println("After sort")
    for ig in 2:6
        @printf("ig = %4d G = [%10.5f,%10.5f,%10.5f] G2 = %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    for ig in Ng-5:Ng
        @printf("ig = %4d G = [%10.5f,%10.5f,%10.5f] G2 = %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end

    display(RecVecs); println()
    println("Pass here in init_gvec_gamma")

    return
end