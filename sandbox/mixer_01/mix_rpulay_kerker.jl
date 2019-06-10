function mix_rpulay_kerker!( pw::PWGrid, x, gx, beta, X, F, iter, MIXDIM, x_old, f_old )

    f = precKerker(pw, gx - x)

    if iter > 1
        dx = x - x_old
        df = f - f_old
    end

# example MIXDIM = 4

# iter = 1, active_dim = 0 , start iteration
# iter = 2, active_dim = 1
# iter = 3, active_dim = 2
# iter = 4, active_dim = 3
# iter = 5, active_dim = 4
# iter = 6, active_dim = 0, restart

    active_dim = mod(iter-1,MIXDIM+1)

    if iter > 1
        # Restart, use the last history
        if active_dim == 0
            #
            X[:,1] = X[:,MIXDIM]
            F[:,1] = F[:,MIXDIM]
            X[:,2:MIXDIM] .= 0.0
            F[:,2:MIXDIM] .= 0.0
            #
            Xk = X[:,1]
            Fk = F[:,1]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            xnew = x + beta*f - addv
        else
            for i = active_dim:-1:2
                X[:,i] = X[:,i-1]
                F[:,i] = F[:,i-1]
            end
            X[:,1] = dx
            F[:,1] = df
            #
            Xk = X[:,1:active_dim]
            Fk = F[:,1:active_dim]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            xnew = x + beta*f - addv
        end
    else
        xnew = x + beta*f
        X[:,1] = x
        F[:,1] = f
    end

    x_old[:] = x
    f_old[:] = f
    return xnew

end