"""
An implementation of periodic Pulay mixing as described in

  Amartya S. Banerjee, Phanish Suryanarayana, John E. Pask
  "Periodic Pulay method for robust and efficient convergence
  accelaration of self-consistent field iterations"
  Chemical Physics Letters 647 (2016) 31-35.
"""
function mix_ppulay!( x, gx, beta, X, F, iter, MIXDIM::Int64, k::Int64,
                      x_old, f_old )

    f = gx - x

    if iter > 1
        dx = x - x_old
        df = f - f_old
    end

    if iter == 1
        is_pulay = false
    else
        is_pulay = iter%k == 0
    end

    active_dim = mod(iter-1,MIXDIM) + 1

    if iter == 1
        X[:,1] = x
        F[:,1] = f
    end

    # shift history
    for i = active_dim:-1:2
        X[:,i] = X[:,i-1]
        F[:,i] = F[:,i-1]
    end

    # new history, put at the beginning
    if iter > 1
        X[:,1] = dx
        F[:,1] = df
    end

    if is_pulay
        Xk = X[:,1:active_dim]
        Fk = F[:,1:active_dim]
        addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
        xnew = x + beta*f - addv
    else
        xnew = x + beta*f
    end

    x_old[:] = x
    f_old[:] = f
    return xnew

end

