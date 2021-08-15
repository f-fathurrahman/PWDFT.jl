function mix_pulay!( x, gx, beta, X, F, iter, MIXDIM::Int64,
                     x_old, f_old )

    f = gx[:] - x[:]

    if iter == 1
        dx = copy(x)
        df = copy(f)
    else
        dx = x[:] - x_old[:]
        df = f[:] - f_old[:]
    end

    # shift history
    if iter <= MIXDIM
        X[:,iter] = dx[:]
        F[:,iter] = df[:]
    else
        for i in 1:MIXDIM-1
            X[:,i] = X[:,i+1]
            F[:,i] = F[:,i+1]
        end
        # new history, put at the end
        X[:,MIXDIM] = dx[:]
        F[:,MIXDIM] = df[:]
    end

    x_old[:] = x[:]
    f_old[:] = f[:]

    # Pulay mixing begins at iter MIXDIM+1
    if iter > MIXDIM
        addv = (X + beta*F)*inv(F'*F)*(F'*f)
        x[:] = x[:] + beta*f - addv
    else
        x[:] = x[:] + beta*f
    end
    
    return

end

