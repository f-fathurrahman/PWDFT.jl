mutable struct BroydenMixer_G
    betamix::Float64
    mixdim::Int64
    df::Vector{Matrix{ComplexF64}}
    dv::Vector{Matrix{ComplexF64}}
    df_bec::Union{Nothing,Vector{Array{Float64}}}
    dv_bec::Union{Nothing,Vector{Array{Float64}}}
    conv_thr::Float64
    is_converged::Bool
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


function BroydenMixer_G(Rhoe::Matrix{ComplexF64}, betamix; mixdim=8, conv_thr=1e-6)
    df = Vector{Matrix{ComplexF64}}(undef, mixdim)
    dv = Vector{Matrix{ComplexF64}}(undef, mixdim)
    for i in 1:mixdim
        df[i] = zeros(ComplexF64, size(Rhoe))
        dv[i] = zeros(ComplexF64, size(Rhoe))
    end
    # No becsum is given
    df_bec = nothing
    dv_bec = nothing
    conv_thr = conv_thr
    is_converged = false
    return BroydenMixer_G(betamix, mixdim, df, dv, df_bec, dv_bec, conv_thr, is_converged)
end


function BroydenMixer_G(
    Rhoe::Matrix{ComplexF64},
    bec::Array{Float64,3},
    betamix;
    mixdim=8, conv_thr=5e-7
)
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
    conv_thr = conv_thr
    is_converged = false
    return BroydenMixer_G(betamix, mixdim, df, dv, df_bec, dv_bec, conv_thr, is_converged)
end


function do_mix!(
    mixer::BroydenMixer_G,
    pw,
    deltain, deltaout_,
    iterSCF::Int64;
    bec_in=nothing, bec_out=nothing
)   
    is_converged = _do_mix_broyden_G!(
        pw,
        deltain, deltaout_,
        mixer.betamix,  # formal parameter name is alphamix
        iterSCF, mixer.mixdim,
        mixer.df, mixer.dv,
        bec_in=bec_in,
        bec_out_=bec_out,
        df_bec=mixer.df_bec,
        dv_bec=mixer.dv_bec,
        conv_thr=mixer.conv_thr
    )
    mixer.is_converged = is_converged
    return
end


function _do_mix_broyden_G!(
    pw,
    deltain, deltaout_,
    alphamix::Float64,
    iterSCF::Int64, n_iter::Int64,
    df, dv;
    bec_in=nothing, bec_out_=nothing,
    df_bec=nothing, dv_bec=nothing,
    conv_thr=5e-7
)

    # Convert to Tot and magn
    Nspin = size(deltain, 2)
    if Nspin == 2
        _convert_to_rhoe_tot_magn!( deltain )
        _convert_to_rhoe_tot_magn!( deltaout_ )
        #@info "Check total charge (G-space) = $(deltain[1,1]*pw.CellVolume)"
        #@info "Check magn (G-space) = $(deltain[1,2]*pw.CellVolume)"
    end

    # df(ndim,n_iter)
    # dv(ndim,n_iter)

    deltaout = copy(deltaout_)  # do not replace deltaout_
    if !isnothing(bec_out_)
        @assert !isnothing(bec_in)
        bec_out = copy(bec_out_)
    end

    deltain_save = copy(deltain)
    if !isnothing(bec_in)
        bec_in_save = copy(bec_in)
    end

    if pw.using_dual_grid
        Ngf = pw.gvecs.Ng
    else
        Ngf = pw.gvec.Ng
    end
    idx_g2r = pw.gvec.idx_g2r

    # call assign_scf_to_mix_type(rhoin, rhoin_m)
    # rhoin_m is deltain

    # call assign_scf_to_mix_type(input_rhout, rhout_m) 
    # rhout_m is deltaout

    # compute differences
    #call mix_type_AXPY( -1.d0, rhoin_m, rhout_m )  ! rhoout_m <- (-1)*rhoin_m + rhout_m
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        deltaout[ip,:] = deltaout[ip,:] - deltain[ip,:]
    end
    if !isnothing(bec_in)
        bec_out[:] = bec_out[:] - bec_in[:]
    end

    dr2 = rhoe_ddot(pw, deltaout, deltaout)
    #@info "iter = $(iterSCF) dr2 = $(dr2)"
    if dr2 < 0.0
        error("dr2 is negative in _do_mix_broyden_G")
    end

    is_converged = false
    if sqrt(dr2) < conv_thr
        @info "sqrt(dr2) = $(sqrt(dr2)) is quite small already. Should be converged"
        is_converged = true
        # Do mixing one more time even if it is already converged
    end

    mixrho_iter = iterSCF # alias
    #
    # iter_used = mixrho_iter-1  if  mixrho_iter <= n_iter
    # iter_used = n_iter         if  mixrho_iter >  n_iter
    #
    iter_used = min( (mixrho_iter - 1), n_iter )
    #
    # ipos is the position in which results from the present iteraction
    # are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    #
    ipos = mixrho_iter - 1 - floor(Int64, (mixrho_iter-2)/n_iter)*n_iter

    inext = mixrho_iter - floor(Int64, (mixrho_iter - 1)/n_iter)*n_iter

    
    # Set up df
    if mixrho_iter > 1
        #@info "ipos = $(ipos) sum df ipos = $(sum(df[ipos][:,1]))"
        #(Nspin == 2) && @info "ipos = $(ipos) sum df ipos = $(sum(df[ipos][:,2]))"
        for ig in 1:Ngf
            ip = idx_g2r[ig]
            #call mix_type_AXPY( -1.d0, rhout_m, df(ipos) )
            df[ipos][ip,:] = (-1)*deltaout[ip,:] + df[ipos][ip,:]
            #call mix_type_AXPY( -1.d0, rhoin_m, dv(ipos) )
            dv[ipos][ip,:] = (-1)*deltain[ip,:] + dv[ipos][ip,:] 
        end
        # also do axpy operation for bec
        if !isnothing(bec_in)
            df_bec[ipos][:] = (-1)*bec_out[:] + df_bec[ipos][:]
            dv_bec[ipos][:] = (-1)*bec_in[:] + dv_bec[ipos][:]  
        end
    end

    # for next iter, save
    deltaout_save = copy(deltaout)
    if !isnothing(bec_in) # check formal arg name
        @assert !isnothing(bec_out_)
        bec_out_save = copy(bec_out)
    end

    if iter_used > 0
    
        beta = zeros(Float64, iter_used, iter_used)
        
        for i in 1:iter_used, j in i:iter_used
            beta[i,j] = rhoe_ddot( pw, df[j], df[i] )
            beta[j,i] = beta[i,j]
        end

        #println("\nbeta matrix before inverse")
        #display(beta[1:iter_used,1:iter_used]); println()    
        beta_inv = inv(beta[1:iter_used,1:iter_used])
        @views beta[1:iter_used,1:iter_used] = beta_inv[:,:]

        work = zeros(Float64, iter_used)
        for i in 1:iter_used
            work[i] = rhoe_ddot( pw, df[i], deltaout )
        end
        
        for i in 1:iter_used
            gamma0 = dot( beta[1:iter_used,i], work[1:iter_used] )
            #@info "gamma0 = $(gamma0)"
            for ig in 1:Ngf
                ip = idx_g2r[ig]
                #call mix_type_AXPY( -gamma0, dv(i), rhoin_m )
                deltain[ip,:] = -gamma0*dv[i][ip,:] + deltain[ip,:]
                #call mix_type_AXPY( -gamma0, df(i), rhout_m )
                deltaout[ip,:] = -gamma0*df[i][ip,:] + deltaout[ip,:] 
            end
            if !isnothing(bec_in)
                bec_in[:] = -gamma0*dv_bec[i][:] + bec_in[:] 
                bec_out[:] = -gamma0*df_bec[i][:] + bec_out[:] 
            end
        end
    end

    #
    # set new trial density
    #
    # rhoin_m <- alphamix*rhout_m + rhoin_m
    #call mix_type_AXPY( alphamix, rhout_m, rhoin_m )
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        deltain[ip,:] =  alphamix*deltaout[ip,:] + deltain[ip,:]
    end
    if !isnothing(bec_in)
        bec_in[:] = alphamix*bec_out[:] + bec_in[:]
    end

    # High freq mixing: using linear mixing
    if pw.using_dual_grid
        ig_start = pw.gvecs.Ng + 1
        ig_stop = pw.gvec.Ng
        for ig in ig_start:ig_stop
            ip = idx_g2r[ig]
            # use original variables (before modification)
            deltain[ip,:] = deltain_save[ip,:] + alphamix *( deltaout_[ip,:] - deltain_save[ip,:] )
        end
    end

    #@info "mixrho_iter = $(mixrho_iter) iter_used = $(iter_used)"
    #@info "ipos = $(ipos) inext = $(inext)"

    # Save data for next iteration
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        df[inext][ip,:] = deltaout_save[ip,:]
        dv[inext][ip,:] = deltain_save[ip,:]
    end
    if !isnothing(bec_in)
        df_bec[inext][:] = bec_out_save[:]
        dv_bec[inext][:] = bec_in_save[:]
    end

    # Convert back to up,dn
    if Nspin == 2
        _convert_to_rhoe_up_dn!( deltain )
        _convert_to_rhoe_up_dn!( deltaout_ ) # not really needed
    end

    return is_converged

end