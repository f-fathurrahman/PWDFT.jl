"""
An implementation of periodic Pulay mixing as described in

  Amartya S. Banerjee, Phanish Suryanarayana, John E. Pask
  "Periodic Pulay method for robust and efficient convergence
  accelaration of self-consistent field iterations"
  Chemical Physics Letters 647 (2016) 31-35.
"""
function mix_ppulay!( x, gx, beta, X, F, iter, MIXDIM::Int64, k::Int64,
                      x_old, f_old )

    f = gx[:] - x[:]

    if iter == 1
        dx = copy(x)
        df = copy(f)
    else
        dx = x[:] - x_old[:]
        df = f[:] - f_old[:]
    end

    if iter < MIXDIM
        is_pulay = false
    else
        is_pulay = iter%k == 0
    end

    # shift history
    if iter <= MIXDIM
        X[:,iter] = dx[:]
        F[:,iter] = df[:]
    else
        for i = 1:MIXDIM-1
            X[:,i] = X[:,i+1]
            F[:,i] = F[:,i+1]
        end
    end

    # new history, put at the end
    if iter > MIXDIM
        X[:,MIXDIM] = dx[:]
        F[:,MIXDIM] = df[:]
    end

    x_old[:] = x[:]
    f_old[:] = f[:]

    if is_pulay
        addv = (X + beta*F)*inv(F'*F)*(F'*f) 
        x[:] = x[:] + beta*f - addv
    else
        x[:] = x[:] + beta*f
    end

    return

end

