#
# Adapted from andersonmix of KSSOLV
# Spin-polarized version
function andersonmix!( vin::Array{Float64,2}, vout::Array{Float64,2},
    beta::Float64, df::Array{Float64,3}, dv::Array{Float64,3},
    iter::Int64, mixdim::Int64)
    # output [vnew,df,dv]

    Nspin = size(vin,2)
    Npoints = size(vin)[1]

    dvout = zeros(Float64,Npoints,Nspin)

    # function evaluation overwrites vout,
    # use dvout instead
    dvout[:,:] = vout[:,:] - vin[:,:]

    iterused = min(iter-1,mixdim)
    ipos = iter - 1 - floor(Int64, (iter-2)/mixdim)*mixdim

    #@printf("ipos = %d\n", ipos)

    if iter > 1
        # compute the changes in function evaluations and the step (changes in potentials)
        df[:,ipos,:] = df[:,ipos,:] - dvout[:,:]
        dv[:,ipos,:] = dv[:,ipos,:] - vin[:,:]
    end

    vinsave = zeros(Float64,Npoints,Nspin)
    dvoutsave = zeros(Float64,Npoints,Nspin)

    vinsave[:,:]  = vin[:,:]
    dvoutsave[:,:] = dvout[:,:]

    if (iter > 1)
        gammas = pinv(df[:,1:iterused])*dvout
        for i = 1:iterused
            vin[:,:]  = vin[:,:]  - gammas[i] * dv[:,i,:]
            dvout[:,:] = dvout[:,:] - gammas[i] * df[:,i,:]
        end
    end

    inext = iter - floor( Int64, (iter - 1) / mixdim) * mixdim
    #@printf("inext = %d\n", inext)

    df[:,inext,:] = dvoutsave[:,:]
    dv[:,inext,:] = vinsave[:,:]

    return vin[:,:] + beta*dvout[:,:]
end



#
# Adapted from andersonmix of KSSOLV
#
function andersonmix!( vin::Array{Float64,1}, vout::Array{Float64,1},
                       beta::Float64, df::Array{Float64,2}, dv::Array{Float64,2},
                       iter::Int64, mixdim::Int64)
    # output [vnew,df,dv]
    Npoints = size(vin)[1]

    dvout = zeros(Float64,Npoints)
    
    # function evaluation overwrites vout,
    # use dvout instead
    dvout[:] = vout[:] - vin[:]

    println("in anderson, beta = ", beta)
    println("sum(vin) = ", sum(vin))
    println("sum(vout) = ", sum(vout))
    println("sum(dvout) = ", sum(dvout))    

    iterused = min(iter-1,mixdim)
    ipos = iter - 1 - floor(Int64, (iter-2)/mixdim)*mixdim

    #@printf("ipos = %d\n", ipos)

    if iter > 1
        # compute the changes in function evaluations and the step (changes in potentials)
        df[:,ipos] = df[:,ipos] - dvout
        dv[:,ipos] = dv[:,ipos] - vin
    end
    
    vinsave = zeros(Float64,Npoints)
    dvoutsave = zeros(Float64,Npoints)

    vinsave[:]  = vin[:]
    dvoutsave[:] = dvout[:]

    if (iter > 1)
        gammas = pinv(df[:,1:iterused])*dvout
        for i = 1:iterused
            vin[:]  = vin[:]  - gammas[i] * dv[:,i]
            dvout[:] = dvout[:] - gammas[i] * df[:,i]
        end
    end

    inext = iter - floor( Int64, (iter - 1) / mixdim) * mixdim
    #@printf("inext = %d\n", inext)

    df[:,inext] = dvoutsave[:]
    dv[:,inext] = vinsave[:]

    return vin[:] + beta*dvout[:]
end
