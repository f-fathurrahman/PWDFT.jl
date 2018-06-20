function rgen( dtau, rmax, mxr, at, bg )
    #
    #   generates neighbours shells (cartesian, in units of lattice parameter)
    #   with length < rmax,and returns them in order of increasing length:
    #      r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
    #   where a1, a2, a3 are primitive lattice vectors. Other input variables:
    #     mxr = maximum number of vectors
    #     at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    #     bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    #   Other output variables:
    #     nrm = the number of vectors with r^2 < rmax^2
    #

    nrm = 0
    r = zeros(Float64,3,mxr)
    r2 = zeros(Float64,mxr)

#  INTEGER, ALLOCATABLE :: irr (:)
#  INTEGER ::  nm1, nm2, nm3, i, j, k, ipol, ir, indsw, iswap
#  real(DP) :: tt, swap
#  real(DP), EXTERNAL :: dnrm2

    ds = zeros(Float64,3)
    dtau0 = zeros(Float64,3)
    t = zeros(Float64,3)

    nrm = 0
    if rmax == 0.0
        return
    end

# bring dtau into the unit cell centered on the origin - prevents trouble
# if atomic positions are not centered around the origin but displaced
# far away (remember that translational invariance allows this!)

    ds[:] = bg[:,:]*dtau[:]
    ds[:] = ds[:] - round.( ds[:] .+ 0.5 )
    dtau0[:] = at[:,:]*ds[:]

    irr = zeros(Int64,mxr)
  
    # these are estimates of the maximum values of needed integer indices

    nm1 = round( Int64, norm(bg[:,1])*rmax ) + 4
    nm2 = round( Int64, norm(bg[:,2])*rmax ) + 4
    nm3 = round( Int64, norm(bg[:,3])*rmax ) + 4

    @printf("rmax = %18.10f\n", rmax)
    @printf("nm = [%d,%d,%d]\n", nm1, nm2, nm3)
  
    for i = -nm1:nm1
    for j = -nm2:nm2
    for k = -nm3:nm3
        tt = 0.0
        for ipol = 1:3
            t[ipol] = i*at[ipol,1] + j*at[ipol,2] + k*at[ipol,3] - dtau0[ipol]
            tt = tt + t[ipol]*t[ipol]
        end
        if (tt <= rmax^2) & (abs(tt) > 1.e-10)
            nrm = nrm + 1
            if nrm > mxr
                @printf("ERROR in rgen: too many r-vectors: %d", nrm)
                exit()
            end
            for ipol = 1:3
                r[ipol,nrm] = t[ipol]
            end
            r2[nrm] = tt
        end
    end
    end
    end
  
    # reorder the vectors in order of increasing magnitude
    if nrm > 1
        irr = sortperm(r2[1:nrm])
        r2[1:nrm] = r2[irr]
        r2[nrm+1:end] .= 0.0
        r[:,1:nrm] = r[1:3,irr]
        r[:,nrm+1:end] .= 0.0
    end

"""
    for ir = 1:nrm - 1
        indsw = irr[ir]
        if indsw /= ir
            for ipol = 1, 3
                swap = r[ipol,indsw]
                r[ipol,indsw] = r[ipol,irr[indsw]]
                r[ipol,irr[indsw]] = swap
            end
            iswap = irr[ir]
            irr[ir] = irr[indsw]
            irr[indsw] = iswap
        end
    end
"""
    @printf("nrm = %d\n", nrm)
    return nrm, r, r2
end

