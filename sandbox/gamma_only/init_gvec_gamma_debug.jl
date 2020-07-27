function inv_mm_to_nn(nn::Int64, S::Int64)
    if nn < 0
        return nn + S
    else
        return nn
    end
    #if mm > S/2
    #    return mm - S
    #else
    #    return mm
    #end
end


#
# Based on ggen of QE-6.5
#

function init_gvec_gamma( Ns, RecVecs, ecutrho )
    
    ni = floor( Int64, (Ns[1]-1)/2 )
    nj = floor( Int64, (Ns[2]-1)/2 )
    nk = floor( Int64, (Ns[3]-1)/2 )
    #
    # gamma-only: exclude space with x < 0
    #
    istart = 0

    tx = zeros(3)
    ty = zeros(3)
    t = zeros(3)
    tt = zeros(Ns[3])

    ngm = 0
    ip = 0
    ipm = 0

    # counter
    ipc = 0
    ipmc = 0

    SMALL = eps()

    for i in istart:ni
       #
       # gamma-only: exclude plane with x = 0, y < 0
       #
       if i == 0
          jstart = 0
       else
          jstart = -nj
       end
       
       tx[1:3] = i * RecVecs[1:3,1] # FIXME: Check this
       
       for j in jstart:nj
          
          #IF ( .NOT. global_sort ) THEN
          #   IF ( fft_stick_index( dfftp, i, j ) == 0 ) CYCLE jloop
          #   is_local = .TRUE.
          #ELSE
          #   IF ( dfftp%lpara .AND. fft_stick_index( dfftp, i, j ) == 0) THEN
          #      is_local = .FALSE.
          #   ELSE
          #      is_local = .TRUE.
          #   END IF
          #END IF
          #
          # gamma-only: exclude line with x = 0, y = 0, z < 0
          #
          if ( (i == 0) && (j == 0) )
             kstart = 0
          else
             kstart = -nk
          end
          
          ty[1:3] = tx[1:3] + j * RecVecs[1:3,2]
          #
          #  compute all the norm square
          #
          for k in kstart:nk
             #
             ip = ip + 1
             t[1] = ty[1] + k * RecVecs[1,3]
             t[2] = ty[2] + k * RecVecs[2,3]
             t[3] = ty[3] + k * RecVecs[3,3]
             tt[k-kstart+1] = t[1]^2 + t[2]^2 + t[3]^2
          end
          #
          #  save all the norm square within cutoff
          #
          for k in kstart:nk
            if 0.5*tt[k-kstart+1] <= ecutrho
                ngm = ngm + 1
                @printf("idx_miller = [%4d %4d %4d] ", i, j, k)
                ip1 = inv_mm_to_nn(i, Ns[1])
                ip2 = inv_mm_to_nn(j, Ns[2])
                ip3 = inv_mm_to_nn(k, Ns[3])
                @printf(" [%4d %4d %4d] ", ip1, ip2, ip3)
                ig1 = PWDFT.mm_to_nn(ip1, Ns[1])
                ig2 = PWDFT.mm_to_nn(ip2, Ns[2])
                ig3 = PWDFT.mm_to_nn(ip3, Ns[3])
                @printf(" [%4d %4d %4d]\n", ig1, ig2, ig3)
                #
                ip  = ip1 + 1 + ip2*Ns[2] + ip3*Ns[2]*Ns[3]
                ipc = ipc + 1
                #
                if ( tt[k-kstart+1] > SMALL )
                    ip1m = inv_mm_to_nn(-i, Ns[1])
                    ip2m = inv_mm_to_nn(-j, Ns[2])
                    ip3m = inv_mm_to_nn(-k, Ns[3])
                    ipm = ip1m + 1 + ip2m*Ns[2] + ip3m*Ns[2]*Ns[3]
                    @printf("Negative   : [%4d %4d %4d]\n", ip1m, ip2m, ip3m)
                    ipmc = ipmc + 1
                end
            end # if
          end # kstart:nk
       end
    end

    println("Ns = ", Ns)
    println("prod(Ns) = ", prod(Ns))
    println("prod(Ns)/2 = ", prod(Ns)/2)
    println("ni*nj*nk = ", ni*nj*nk)
    @printf("ni,nj,nk =  %d, %d, %d\n", ni, nj, nk)
    println("ip last  = ", ip)
    println("ipm last = ", ipm)
    println("ngm = ", ngm)
    println("ipc  = ", ipc)
    println("ipmc = ", ipmc)

    return
end