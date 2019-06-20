#
# Adapted from PWSCF/Yambo
#
function mix_broyden!( deltain, deltaout, alphamix::Float64, iter::Int64, n_iter::Int64, df, dv )
    # df(ndim,n_iter)
    # dv(ndim,n_iter)

    #deltaout = copy(deltaout_)  # do not replace deltaout_

    maxter = 8
    wg0 = 0.01
    wg = ones(maxter)

    deltainsave = copy( deltain )
    #
    # iter_used = iter-1  IF iter <= n_iter
    # iter_used = n_iter  IF iter >  n_iter
    #
    iter_used = min(iter-1,n_iter)
    #
    # ipos is the position in which results from the present iteraction
    # are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    #
    ipos = iter - 1 - Int64(floor((iter-2)/n_iter))*n_iter

    deltaout[:] = deltaout[:] - deltain[:]

    if iter > 1
        df[:,ipos] = deltaout[:] - df[:,ipos]
        dv[:,ipos] = deltain[:]  - dv[:,ipos]
        nrm = sqrt( norm(df[:,ipos])^2 )
        df[:,ipos] = df[:,ipos]/nrm
        dv[:,ipos] = dv[:,ipos]/nrm
    end

    beta = zeros(maxter,maxter)

    for i = 1:iter_used
        for j = i+1:iter_used
            beta[i,j] = wg[i] * wg[j] * dot( df[:,j], df[:,i] )
            beta[j,i] = beta[i,j]
        end
        beta[i,i] = wg0^2 + wg[i]^2
    end

    beta_inv = inv(beta[1:iter_used,1:iter_used])

    beta[1:iter_used,1:iter_used] = beta_inv[:,:]

    work = zeros(iter_used)
    for i = 1:iter_used
        work[i] = dot(df[:,i], deltaout)
    end
    
    deltain[:] = deltain[:] + alphamix*deltaout[:]

    for i = 1:iter_used
        gammamix = 0.0
        for j = 1:iter_used
            gammamix = gammamix + beta[j,i] * wg[j] * work[j]
        end
        deltain[:] = deltain[:] - wg[i]*gammamix*( alphamix*df[:,i] + dv[:,i] )
    end

    inext = iter - Int64(floor((iter - 1)/n_iter))*n_iter
    df[:,inext] = deltaout[:]
    dv[:,inext] = deltainsave[:]

    return

end