mutable struct BroydenMixer_G
    betamix::Float64
    mixdim::Int64
    df::Vector{Matrix{ComplexF64}}
    dv::Vector{Matrix{ComplexF64}}
    df_bec::Union{Nothing,Vector{Array{Float64}}}
    dv_bec::Union{Nothing,Vector{Array{Float64}}}
end

function _convert_to_rhoe_tot_magn!( Rhoe )
    Nspin = size(Rhoe, 2)
    @assert Nspin == 2
    Npoints = size(Rhoe, 1)
    Rhoe_tmp = copy(Rhoe)
    for ip in 1:Npoints
        Rhoe[ip,1] = Rhoe_tmp[ip,1] + Rhoe_tmp[ip,2]
        Rhoe[ip,2] = Rhoe_tmp[ip,1] - Rhoe_tmp[ip,2] 
    end
    return
end

#=
Convert to Rhoe_up and Rhoedn    

Rhoe_tot = Rhoe_up + Rhoe_dn   (A)
magn     = Rhoe_up - Rhoe_dn   (B)

Sum: (A) + (B)
Rhoe_tot + magn = 2*Rhoe_up -> Rhoe_up = 0.5*( Rhoe_tot + magn )

Diff: (A) - (B)
Rhoe_tot - magn = 2*Rhoe_dn -> Rhoe_dn = 0.5*( Rhoe_tot - magn )
=#

function _convert_to_rhoe_up_dn!( Rhoe )
    Nspin = size(Rhoe, 2)
    @assert Nspin == 2
    Npoints = size(Rhoe, 1)
    Rhoe_tmp = copy(Rhoe)
    for ip in 1:Npoints
        Rhoe[ip,1] = 0.5*(Rhoe_tmp[ip,1] + Rhoe_tmp[ip,2])
        Rhoe[ip,2] = 0.5*(Rhoe_tmp[ip,1] - Rhoe_tmp[ip,2])
    end
    return
end


function BroydenMixer_G(Rhoe::Matrix{ComplexF64}, betamix; mixdim=8)
    df = Vector{Matrix{ComplexF64}}(undef, mixdim)
    dv = Vector{Matrix{ComplexF64}}(undef, mixdim)
    for i in 1:mixdim
        df[i] = zeros(ComplexF64, size(Rhoe))
        dv[i] = zeros(ComplexF64, size(Rhoe))
    end
    # No becsum is given
    df_bec = nothing
    dv_bec = nothing
    return BroydenMixer_G(betamix, mixdim, df, dv, df_bec, dv_bec)
end


function BroydenMixer_G(Rhoe::Matrix{ComplexF64}, bec::Array{Float64,3}, betamix; mixdim=8)
    df = Vector{Matrix{ComplexF64}}(undef, mixdim)
    dv = Vector{Matrix{ComplexF64}}(undef, mixdim)
    for i in 1:mixdim
        df[i] = zeros(ComplexF64, size(Rhoe))
        dv[i] = zeros(ComplexF64, size(Rhoe))
    end
    df_bec = Vector{Array{Float64,3}}(undef, mixdim)
    dv_bec = Vector{Array{Float64,3}}(undef, mixdim)
    for i in 1:mixdim
        df_bec[i] = zeros(Float64, size(bec))
        dv_bec[i] = zeros(Float64, size(bec))
    end
    return BroydenMixer_G(betamix, mixdim, df, dv, df_bec, dv_bec)
end


function do_mix!(
    mixer::BroydenMixer_G,
    pw,
    deltain, deltaout_,
    iterSCF::Int64;
    bec_in=nothing, bec_out=nothing
)
    # XXX: alphamix -> mixer.betamix for now
    _do_mix_broyden_G!(
        pw,
        deltain, deltaout_,
        mixer.betamix,
        iterSCF, mixer.mixdim,
        mixer.df, mixer.dv,
        bec_in=bec_in,
        bec_out_=bec_out,
        df_bec=mixer.df_bec,
        dv_bec=mixer.dv_bec
    )
    return
end


#
# Adapted from PWSCF/Yambo
#
# XXX: rename variable alphamix to betamix?
# OUTPUT: deltain
# deltaout_ will not be modified
function _do_mix_broyden_G!(
    pw,
    deltain, deltaout_,
    alphamix::Float64,
    iter::Int64, n_iter::Int64,
    df, dv;
    bec_in=nothing, bec_out_=nothing,
    df_bec=nothing, dv_bec=nothing
)

    # Convert to Tot and magn
    Nspin = size(deltain, 2)
    if Nspin == 2
        _convert_to_rhoe_tot_magn!( deltain )
        _convert_to_rhoe_tot_magn!( deltaout_ )
        @info "Check total charge (G-space) = $(deltain[1,1]*pw.CellVolume)"
        @info "Check magn (G-space) = $(deltain[1,2]*pw.CellVolume)"
    end

    # df(ndim,n_iter)
    # dv(ndim,n_iter)

    deltaout = copy(deltaout_)  # do not replace deltaout_
    if !isnothing(bec_out_)
        @assert !isnothing(bec_in)
        bec_out = copy(bec_out_)
    end

    maxter = n_iter  # FIXME: should be checked against calling functions

    deltainsave = copy(deltain)
    if !isnothing(bec_in)
        bec_in_save = copy(bec_in)
    end
    #
    # iter_used = iter-1  IF iter <= n_iter
    # iter_used = n_iter  IF iter >  n_iter
    #
    iter_used = min(iter-1, n_iter)
    #
    # ipos is the position in which results from the present iteraction
    # are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    #
    ipos = iter - 1 - floor(Int64, (iter-2)/n_iter)*n_iter

    #println("mix_broyden: ipos = ", ipos)
    if pw.using_dual_grid
        Ngf = pw.gvecs.Ng
    else
        Ngf = pw.gvec.Ng
    end
    idx_g2r = pw.gvec.idx_g2r

    for ig in 1:Ngf
        ip = idx_g2r[ig]
        deltaout[ip,:] = deltaout[ip,:] - deltain[ip,:]
    end
    @printf("mix rhoe = %4d, norm diff deltaout[:,1] = %15.10e\n", iter, norm(deltaout[:,1]))
    (Nspin == 2) && @printf("mix rhoe = %4d, norm diff deltaout[:,2] = %15.10e\n", iter, norm(deltaout[:,2]))

    if !isnothing(bec_in)
        bec_out[:] = bec_out[:] - bec_in[:]
        @printf("mix rhoe = %4d, norm diff bec_out[:,:,1] = %15.10e\n", iter, norm(bec_out[:,:,1]))
        (Nspin == 2) && @printf("mix rhoe = %4d, norm diff bec_out[:,:,2] = %15.10e\n", iter, norm(bec_out[:,:,2]))    
    end

    if iter > 1
        for ig in 1:Ngf
            ip = idx_g2r[ig]
            df[ipos][ip,:] = deltaout[ip,:] - df[ipos][ip,:]
            dv[ipos][ip,:] = deltain[ip,:]  - dv[ipos][ip,:]
        end
        nrm = sqrt(rhoe_ddot(pw, df[ipos], df[ipos]))
        @views df[ipos][:,:] = df[ipos][:,:]/nrm
        @views dv[ipos][:,:] = dv[ipos][:,:]/nrm
        #
        if !isnothing(bec_in)
            df_bec[ipos][:] = bec_out[:] - df_bec[ipos][:]
            dv_bec[ipos][:] = bec_in[:]  - dv_bec[ipos][:]
            nrm = norm(df_bec[ipos])
            @views df_bec[ipos][:] = df_bec[ipos][:]/nrm
            @views dv_bec[ipos][:] = dv_bec[ipos][:]/nrm
        end
        # No need to normalize ?
    end

    beta = zeros(Float64, maxter, maxter)

    for i in 1:iter_used
        for j in i+1:iter_used
            beta[i,j] = rhoe_ddot(pw, df[j], df[i])
            beta[j,i] = beta[i,j]
        end
        beta[i,i] = 1.01
    end

    println("\nbeta matrix before inverse")
    display(beta[1:iter_used,1:iter_used]); println()

    beta_inv = inv(beta[1:iter_used,1:iter_used])
    @views beta[1:iter_used,1:iter_used] = beta_inv[:,:]
    
    println("\nbeta matrix after inverse")
    display(beta[1:iter_used,1:iter_used]); println()

    work = zeros(Float64, iter_used)
    for i in 1:iter_used
        work[i] = rhoe_ddot(pw, df[i], deltaout)
    end

    for ig in 1:Ngf
        ip = idx_g2r[ig]
        deltain[ip,:] = deltain[ip,:] + alphamix*deltaout[ip,:]
    end
    if !isnothing(bec_in)
        bec_in[:] = bec_in[:] + alphamix*bec_out[:]
    end

    for i in 1:iter_used
        gammamix = 0.0
        for j in 1:iter_used
            gammamix = gammamix + beta[j,i] * work[j]
        end
        @info "gammamix = $(gammamix)" # can be negative?
        for ig in 1:Ngf
            ip = idx_g2r[ig]
            deltain[ip,:] = deltain[ip,:] - gammamix*( alphamix*df[i][ip,:] + dv[i][ip,:] )
        end
        if !isnothing(bec_in)
            bec_in[:] = bec_in[:] - gammamix*( alphamix*df_bec[i][:] + dv_bec[i][:] )
        end
    end

    # High freq mixing: using linear mixing
    if pw.using_dual_grid
        ig_start = pw.gvecs.Ng + 1
        ig_stop = pw.gvec.Ng
        for ig in ig_start:ig_stop
            ip = idx_g2r[ig]
            deltain[ip,:] = deltain[ip,:] + alphamix *( deltaout[ip,:] - deltain[ip,:] )
        end
    end

    inext = iter - floor(Int64, (iter - 1)/n_iter)*n_iter
    #println("inext = ", inext)
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        df[inext][ip,:] = deltaout[ip,:]
        dv[inext][ip,:] = deltainsave[ip,:]
    end
    if !isnothing(bec_in)
        df_bec[inext][:] = bec_out[:]
        dv_bec[inext][:] = bec_in_save[:]
    end

    # Convert back to up,dn
    if Nspin == 2
        _convert_to_rhoe_up_dn!( deltain )
        _convert_to_rhoe_up_dn!( deltaout_ )
    end

    return

end