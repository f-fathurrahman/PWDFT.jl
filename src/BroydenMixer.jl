mutable struct BroydenMixer
    betamix::Float64
    mixdim::Int64
    df::Matrix{Float64}
    dv::Matrix{Float64}
end


function BroydenMixer(Rhoe::Array{Float64}, betamix; mixdim=8)
    df = zeros(Float64, length(Rhoe), mixdim)
    dv = zeros(Float64, length(Rhoe), mixdim)
    return BroydenMixer(betamix, mixdim, df, dv)
end


function do_mix!(
    mixer::BroydenMixer,
    deltain, deltaout_,
    iterSCF::Int64
)
    # XXX: alphamix -> mixer.betamix for now
    _do_mix_broyden!(
        deltain, deltaout_,
        mixer.betamix,
        iterSCF, mixer.mixdim,
        mixer.df, mixer.dv
    )
    return
end


#
# Adapted from PWSCF/Yambo
#
# XXX: rename variable alphamix to betamix?
# OUTPUT: deltain
# deltaout_ will not be modified
function _do_mix_broyden!(
    deltain, deltaout_,
    alphamix::Float64,
    iter::Int64, n_iter::Int64,
    df, dv 
)
    # df(ndim,n_iter)
    # dv(ndim,n_iter)

    deltaout = copy(deltaout_)  # do not replace deltaout_

    maxter = n_iter  # FIXME: should be checked against calling functions
    wg0 = 0.01
    #wg0 = 0.0 # QE
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
    ipos = iter - 1 - floor(Int64, (iter-2)/n_iter)*n_iter

    #println("mix_broyden: ipos = ", ipos)

    @views deltaout[:] = deltaout[:] - deltain[:]
    #@printf("iter = %4d, norm diff Rhoe (deltaout) = %15.10e\n", iter, norm(deltaout))

    if iter > 1
        @views df[:,ipos] = deltaout[:] - df[:,ipos]
        @views dv[:,ipos] = deltain[:]  - dv[:,ipos]
        @views nrm = norm(df[:,ipos])
        @views df[:,ipos] = df[:,ipos]/nrm
        @views dv[:,ipos] = dv[:,ipos]/nrm
    end

    beta = zeros(maxter,maxter)

    for i in 1:iter_used
        for j in i+1:iter_used
            beta[i,j] = wg[i] * wg[j] * real(dot(df[:,j],df[:,i]))
            beta[j,i] = beta[i,j]
        end
        beta[i,i] = wg0^2 + wg[i]^2
    end

    #println("\nbeta matrix before inverse")
    #display(beta); println()

    beta_inv = inv(beta[1:iter_used,1:iter_used])

    @views beta[1:iter_used,1:iter_used] = beta_inv[:,:]
    
    #println("\nbeta matrix after inverse")
    #display(beta[1:iter_used,1:iter_used]); println()

    work = zeros(iter_used)
    for i in 1:iter_used
        work[i] = real(dot(df[:,i], deltaout))
    end
    
    @views deltain[:] = deltain[:] + alphamix*deltaout[:]

    for i in 1:iter_used
        gammamix = 0.0
        for j in 1:iter_used
            gammamix = gammamix + beta[j,i] * wg[j] * work[j]
        end
        @views deltain[:] = deltain[:] - wg[i]*gammamix*( alphamix*df[:,i] + dv[:,i] )
    end

    inext = iter - floor(Int64, (iter - 1)/n_iter)*n_iter
    #println("inext = ", inext)

    @views df[:,inext] = deltaout[:]
    @views dv[:,inext] = deltainsave[:]

    return

end

function do_mix_precKerker!(
    mixer::BroydenMixer,
    pw::PWGrid,
    deltain, deltaout_,
    iterSCF::Int64
)
    _do_mix_broyden!(
        pw,
        deltain, deltaout_,
        mixer.betamix,
        iterSCF, mixer.mixdim,
        mixer.df, mixer.dv
    )
    return
end


function _do_mix_broyden!(
    pw::PWGrid,
    deltain, deltaout_,
    alphamix::Float64,
    iter::Int64, n_iter::Int64,
    df, dv 
)
    deltaout = copy(deltaout_)  # do not replace deltaout_
    maxter = n_iter
    wg0 = 0.01
    wg = ones(maxter)

    deltainsave = copy( deltain )
    iter_used = min(iter-1,n_iter)
    ipos = iter - 1 - floor(Int64, (iter-2)/n_iter)*n_iter

    deltaout[:] = precKerker(pw, deltaout - deltain)

    if iter > 1
        @views df[:,ipos] = deltaout[:] - df[:,ipos]
        @views dv[:,ipos] = deltain[:]  - dv[:,ipos]
        @views nrm = norm(df[:,ipos])
        @views df[:,ipos] = df[:,ipos]/nrm
        @views dv[:,ipos] = dv[:,ipos]/nrm
    end

    beta = zeros(maxter,maxter)

    for i in 1:iter_used
        for j in i+1:iter_used
            beta[i,j] = wg[i] * wg[j] * real(dot(df[:,j],df[:,i]))
            beta[j,i] = beta[i,j]
        end
        beta[i,i] = wg0^2 + wg[i]^2
    end

    beta_inv = inv(beta[1:iter_used,1:iter_used])
    work = zeros(iter_used)
    for i in 1:iter_used
        work[i] = real(dot(df[:,i], deltaout))
    end
    
    @views deltain[:] = deltain[:] + alphamix*deltaout[:]

    for i in 1:iter_used
        gammamix = 0.0
        for j in 1:iter_used
            gammamix = gammamix + beta[j,i] * wg[j] * work[j]
        end
        @views deltain[:] = deltain[:] - wg[i]*gammamix*( alphamix*df[:,i] + dv[:,i] )
    end

    inext = iter - floor(Int64, (iter - 1)/n_iter)*n_iter
    @views df[:,inext] = deltaout[:]
    @views dv[:,inext] = deltainsave[:]
    return
end
