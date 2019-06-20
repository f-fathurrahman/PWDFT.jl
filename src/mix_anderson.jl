#
# This function is adapted from Anderson mixing function in KSSOLV
#
# vin is the result
# vin will be overwritten for next SCF iteration
function mix_anderson!( vin, vout,
                        beta::Float64, df::Array{Float64,2}, dv::Array{Float64,2},
                        iter::Int64, mixdim::Int64 )
    # Residual
    dvout = vout[:] - vin[:]

    iterused = min(iter-1,mixdim)
    ipos = iter - 1 - floor(Int64, (iter-2)/mixdim)*mixdim

    if iter > 1
        df[:,ipos] = df[:,ipos] - dvout[:]
        dv[:,ipos] = dv[:,ipos] - vin[:]
    end

    vinsave  = copy(vin)
    dvoutsave = copy(dvout)

    if (iter > 1)
        gammas = pinv(df[:,1:iterused])*dvout  
        for i = 1:iterused
            vin[:]  = vin[:]  - gammas[i] * dv[:,i]
            dvout[:] = dvout[:] - gammas[i] * df[:,i]
        end
    end

    inext = iter - floor( Int64, (iter - 1) / mixdim) * mixdim

    df[:,inext] = dvoutsave[:]
    dv[:,inext] = vinsave[:]

    vin[:] = vin[:] + beta*dvout[:]

    return
end

