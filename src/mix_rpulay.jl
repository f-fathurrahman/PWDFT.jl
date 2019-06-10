"""
An implementation of restarted Pulay mixing as described in

  Phanisri P. Pratapa and Phanish Suryanarayana.
  Restarted Pulay mixing for efficient and robust acceleration of
  fixed-point iterations.
  Chemical Physics Letters 635 (2015) 69–74.
"""

# Result (new density) will be written in x
function mix_rpulay!( x, gx, beta, X, F, iter, MIXDIM, x_old, f_old )

    Npts = length(x)
    f = zeros(Float64,Npts)

    for i in 1:Npts
        f[i] = gx[i] - x[i]
    end

    if iter > 1
        dx = zeros(Float64,Npts)
        df = zeros(Float64,Npts)
        for i in 1:Npts
            dx[i] = x[i] - x_old[i]
            df[i] = f[i] - f_old[i]
        end
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
            X[:,1] = @view X[:,MIXDIM]
            F[:,1] = @view F[:,MIXDIM]
            X[:,2:MIXDIM] .= 0.0
            F[:,2:MIXDIM] .= 0.0
            #
            Xk = @view X[:,1]
            Fk = @view F[:,1]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            for i = 1:Npts
                # probably x_old and f_old is not needed
                x_old[i] = x[i]  
                f_old[i] = f[i]
                x[i] = x[i] + beta*f[i] - addv[i]
            end

        else

            for i = active_dim:-1:2
                X[:,i] = @view X[:,i-1]
                F[:,i] = @view F[:,i-1]
            end
            X[:,1] = dx
            F[:,1] = df
            #
            Xk = @view X[:,1:active_dim]
            Fk = @view F[:,1:active_dim]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            for i = 1:Npts
                # probably x_old and f_old is not needed
                x_old[i] = x[i]  
                f_old[i] = f[i]
                x[i] = x[i] + beta*f[i] - addv[i]
            end
        end
    else
        for i in 1:Npts
            X[i,1] = x[i]
            F[i,1] = f[i]
            # probably x_old and f_old is not needed
            x_old[i] = x[i]  
            f_old[i] = f[i]
            #
            x[i] = x[i] + beta*f[i]  # write the result
        end
    end
    return
end
