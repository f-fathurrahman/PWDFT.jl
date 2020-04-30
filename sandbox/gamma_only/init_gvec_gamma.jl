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
                
                #IF (ngm > ngm_max) CALL errore ('ggen 1', 'too many g-vectors', ngm)
                #IF ( tt(k-kstart+1) > eps8 ) THEN
                #   g2sort_g(ngm) = tt(k-kstart+1)
                #ELSE
                #   g2sort_g(ngm) = 0.d0
                #ENDIF
                #IF (is_local) THEN
                #  ngm_local = ngm_local + 1
                #  mill_unsorted( :, ngm_local ) = (/ i,j,k /)
                #  g2l(ngm) = ngm_local
                #ELSE
                #  g2l(ngm) = 0
                #ENDIF
            end # if
          end # kstart:nk
       end
    end

    println("ni = ", ni)
    println("nj = ", nj)
    println("nk = ", nk)
    println("ip = ", ip)
    println("ngm = ", ngm)

end