#
# adapted from rgen.f90 of QE-6.6
#
function gen_neighbor_shells!(
    dtau, rmax::Float64,
    LatVecs, RecVecs,
    r, r2
)
    #
    # generates neighbours shells (cartesian, in units of lattice parameter)
    # with length < rmax,and returns them in order of increasing length:
    #    r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
    # where a1, a2, a3 are primitive lattice vectors. Other input variables:
    #   mxr = maximum number of vectors
    #   at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    #   bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    # Other output variables:
    #   nrm = the number of vectors with r^2 < rmax^2

    mxr = size(r2,1)

    at = LatVecs
    bg = RecVecs
    nrm = 0 # counter of vectors that are searched
    SMALL = eps()
    if rmax <= SMALL
        return nrm
    end
  
    ds = zeros(Float64, 3)
    dtau0 = zeros(Float64, 3)
    t = zeros(Float64, 3)

    # bring dtau into the unit cell centered on the origin - prevents trouble
    # if atomic positions are not centered around the origin but displaced
    # far away (remember that translational invariance allows this!)
    ds[:] = bg' * dtau
    ds[:] = ds[:] .- round.(ds)
    dtau0[:] = at * ds

    # these are estimates of the maximum values of needed integer indices
    nm1 = floor(Int64, norm(bg[:,1]) * rmax / (2*pi) ) + 2
    nm2 = floor(Int64, norm(bg[:,2]) * rmax / (2*pi) ) + 2
    nm3 = floor(Int64, norm(bg[:,3]) * rmax / (2*pi) ) + 2

    #@printf("nm1 = %5d, nm2 = %5d, nm3 = %5d\n", nm1, nm2, nm3)

    for i in -nm1:nm1, j in -nm2:nm2, k in -nm3:nm3
        tt = 0.0
        for ipol in 1:3
            t[ipol] = i*at[ipol,1] + j*at[ipol,2] + k*at[ipol,3] - dtau0[ipol]
            tt = tt + t[ipol]^2
        end
        if (tt <= rmax^2) && (abs(tt) > 1.e-10)
            nrm = nrm + 1
            if nrm > mxr
                error("too many r-vectors, nrm = $(nrm)")
            end
            for ipol in 1:3
                r[ipol,nrm] = t[ipol]
            end
            r2[nrm] = tt
        end
    end
  
    # reorder the vectors in order of increasing magnitude
    if nrm > 1
        @views irr = sortperm(r2[1:nrm])
        @views r[1:3,1:nrm] .= r[1:3,irr]
        @views r2[1:nrm] .= r2[irr]
        # Set the remaning vectors to zeros (they should not be used)
        @views r[1:3,nrm+1:end] .= 0.0
        @views r2[nrm+1:end] .= 0.0
    end
    return nrm
end

