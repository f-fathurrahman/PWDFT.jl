
#
# Adapted from andersonmix of KSSOLV
#
function andersonmix!( vin::Array{Float64,1}, vout::Array{Float64,1},
                       beta::Float64, df, dv, iter::Int64, mixdim::Int64)
    # output [vnew,df,dv]
    Npoints = size(vin)[1]

    # function evaluation overwrites vout,
    # use dvout instead
    dvout = vout - vin

    iterused = min(iter-1,mixdim)
    ipos = iter - 1 - floor(Int64, (iter-2)/mixdim)*mixdim

    #@printf("ipos = %d\n", ipos)

    if iter > 1
        # compute the changes in function evaluations and the step (changes in potentials)
        df[:,ipos] = df[:,ipos] - dvout
        dv[:,ipos] = dv[:,ipos] - vin
    end
    #
    vinsave  = copy(vin)
    dvoutsave = copy(dvout)

    if (iter > 1)
        #B = (df(:,1:iterused))'*(df(:,1:iterused));
        #B = (B+B')/2;
        #
        #work = df(:,1:iterused)'*reshape(vout,n123,1);
        #gammas = B\(work);
        #
        # the use of pseudo-inverse is not the most efficient
        # implementation.
        #
        gammas = pinv(df[:,1:iterused])*dvout
        #println("gammas = ", gammas[:])
        #
        # the following loop can be replaced by two simple gemvs
        # and some reshape operations.
        #
        for i = 1:iterused
            vin  = vin  - gammas[i] * dv[:,i]
            dvout = dvout - gammas[i] * df[:,i]
        end
    end

    inext = iter - floor( Int64, (iter - 1) / mixdim) * mixdim
    #@printf("inext = %d\n", inext)

    df[:,inext] = dvoutsave
    dv[:,inext] = vinsave

    return vin[:] + beta*dvout[:]
end
