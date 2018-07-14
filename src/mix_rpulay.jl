function mix_rpulay!( x, gx, beta, X, F, iter, MIXDIM, dRhoe )

    f = gx - x
    
    m = MIXDIM
    k = iter - 1
    
# iter = 1, active_dim = 1
# iter = 2, active_dim = 2
# ...
# iter = MIXDIM, active_dim = MIXDIM
# iter = MIXDIM+1, active_dim = 1
# iter = MIXDIM+2, active_dim = 2
# ...

    active_dim = mod(iter,MIXDIM)
    if active_dim == 0
        active_dim = 1
    end

    println("active_dim = ", active_dim)

    if iter > 1
        # Restart, use the last history
        if active_dim == 1
            X[:,1] = X[:,MIXDIM]
            F[:,1] = F[:,MIXDIM]
            X[:,:] .= 0.0
            F[:,:] .= 0.0
        else
            X[:,active_dim] = x - X[:,active_dim-1]
            F[:,active_dim] = f - F[:,active_dim-1]
        end
    else
        xnew = x + beta*f
        X[:,1] = xnew - x
        F[:,1] = f
    end

    return xnew

end