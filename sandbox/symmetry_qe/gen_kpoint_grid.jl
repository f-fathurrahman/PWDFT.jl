# Adapted from Quantum ESPRESSO 6.4: PW/src/kpoint_grid.f90
# 
# Copyright (C) 2001-2007 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License.

function gen_kpoint_grid( LatVecs, nk1, nk2, nk3, k1, k2, k3, time_reversal, s )

    SMALL = 1.0e-5

    RecVecs = 2*pi*inv(LatVecs')

    Nrots = size(s,3)

    nkr = nk1*nk2*nk3
    xkg = zeros(Float64,3, nkr)
    wkk = zeros(Float64, nkr)
    
    xkr = zeros(Float64,3)
    equiv = zeros(Int64, nkr)
    

    for i in 1:nk1, j in 1:nk2, k in 1:nk3
        ik = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
        #  xkg are the components of the complete grid in crystal axis
        xkg[1,ik] = (i-1)/nk1 + (k1)/2.0/nk1
        xkg[2,ik] = (j-1)/nk2 + (k2)/2.0/nk2
        xkg[3,ik] = (k-1)/nk3 + (k3)/2.0/nk3
    end

    #  equiv(nk)  = nk : k-point nk is not equivalent to any previous k-point
    #  equiv(nk) != nk : k-point nk is equivalent to k-point equiv(nk)

    for ik in 1:nkr
        equiv[ik] = ik
    end

    for ik in 1:nkr
        
        # check if this k-point has already been found equivalent to another
        if equiv[ik] == ik

            wkk[ik] = 1.0
            
            # check if there are equivalent k-point to this in the list
            # (excepted those previously found to be equivalent to another)
            # check both k and -k
            for irot in 1:Nrots
                
                for i in 1:3
                    xkr[i] = s[i,1,irot]*xkg[1,ik] + s[i,2,irot]*xkg[2,ik] + s[i,3,irot]*xkg[3,ik]
                    xkr[i] = xkr[i] - round( xkr[i] )
                end
                
                xx = xkr[1]*nk1 - 0.5*k1
                yy = xkr[2]*nk2 - 0.5*k2
                zz = xkr[3]*nk3 - 0.5*k3
                
                in_the_list = (abs(xx-round(xx)) <= SMALL) &&
                              (abs(yy-round(yy)) <= SMALL) &&
                              (abs(zz-round(zz)) <= SMALL)
                
                if in_the_list
                    i = round(Int64, xkr[1]*nk1 - 0.5*k1 + 2*nk1) % nk1 + 1
                    j = round(Int64, xkr[2]*nk2 - 0.5*k2 + 2*nk2) % nk2 + 1
                    k = round(Int64, xkr[3]*nk3 - 0.5*k3 + 2*nk3) % nk3 + 1
                    
                    n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                    if (n > ik) && (equiv[n] == n)
                       equiv[n] = ik
                       wkk[ik] = wkk[ik] + 1.0
                    else
                        if (equiv[n] != ik) || (n < ik)
                            error("kpoint_grid something wrong in the checking algorithm")
                        end
                    end
                end

                if time_reversal
                    xx = -xkr[1]*nk1 - 0.5*k1
                    yy = -xkr[2]*nk2 - 0.5*k2
                    zz = -xkr[3]*nk3 - 0.5*k3
                    in_the_list = ( abs(xx - round(xx)) <= SMALL ) &&
                                  ( abs(yy - round(yy)) <= SMALL ) &&
                                  ( abs(zz - round(zz)) <= SMALL )
                    
                    if in_the_list
                        #
                        i = round(Int64, -xkr[1]*nk1 - 0.5*k1 + 2*nk1) % nk1 + 1
                        j = round(Int64, -xkr[2]*nk2 - 0.5*k2 + 2*nk2) % nk2 + 1
                        k = round(Int64, -xkr[3]*nk3 - 0.5*k3 + 2*nk3) % nk3 + 1
                        
                        n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                        
                        if (n > ik) && (equiv[n] == n)
                            equiv[n] = ik
                            wkk[ik] = wkk[ik] + 1.0
                        else
                            if (equiv[n] != ik) || (n < ik)
                                error("kpoint_grid something wrong in the checking algorithm")
                            end
                        end
                    end # if in_the_list

                end # if time_reversal
            
            end

        end # for
    end # if

    # count irreducible points and order them
    nks = 0

    # Calculate nks
    for ik in 1:nkr
        if equiv[ik] == ik
            nks = nks + 1
        end
    end
    println("nks = ", nks)

    wk = zeros(Float64, nks)
    xk = zeros(Float64, 3, nks)

    nks = 0
    fact = 0.0
    for ik in 1:nkr
        if equiv[ik] == ik
            nks = nks + 1
            wk[nks] = wkk[ik]
            fact = fact + wk[nks]
            # bring back into to the first BZ
            for i in 1:3
                xk[i,nks] = xkg[i,ik] - round(xkg[i,ik])
            end
        end
    end

    # go to cartesian axis (in units 2pi/a0)
    #CALL cryst_to_cart(nks,xk,bg,1)
  
    # normalize weights to one
    wk[:] = wk[:]/fact
    xk[:,:] = RecVecs*xk[:,:]

    for ik in 1:nks
        @printf("%18.10f %18.10f %18.10f: %18.10f\n", xk[1,ik], xk[2,ik], xk[3,ik], wk[ik])
    end

    return

end