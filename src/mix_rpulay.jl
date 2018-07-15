function mix_rpulay!( x, gx, beta, X, F, iter, MIXDIM, x_old, f_old )

    f = gx - x

    println("iter = ", iter)
    println("norm f = ", norm(f))

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
# iter = 6, active_dim = 0, restart ?

    active_dim = mod(iter-1,MIXDIM+1)

    println("\nin mix_rpulay:")
    println("active_dim = ", active_dim)

    if iter > 1
        # Restart, use the last history
        if active_dim == 0
            #
            println("Restart rpulay")
            X[:,1] = X[:,MIXDIM]
            F[:,1] = F[:,MIXDIM]
            X[:,2:MIXDIM] .= 0.0
            F[:,2:MIXDIM] .= 0.0
            #
            Xk = X[:,1]
            Fk = F[:,1]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            xnew = x + beta*f - addv
            println("norm addv = ", norm(addv))
            println("norm xnew = ", norm(xnew))
        else
            for i = active_dim:-1:2
                X[:,i] = X[:,i-1]
                F[:,i] = F[:,i-1]
                println(i, "norm X[:,i] = ", norm(X[:,i]))
                println(i, "norm F[:,i] = ", norm(F[:,i]))
            end
            X[:,1] = dx
            F[:,1] = df
            println("norm X[:,1] = ", norm(X[:,1]))
            println("norm F[:,1] = ", norm(F[:,1]))
            #
            Xk = X[:,1:active_dim]
            Fk = F[:,1:active_dim]
            addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
            xnew = x + beta*f - addv
            println("norm addv = ", norm(addv))
            println("norm xnew = ", norm(xnew))
        end
    else
        xnew = x + beta*f
        X[:,1] = x
        F[:,1] = f
        println("norm X[:,1] = ", norm(X[:,1]))
        println("norm F[:,1] = ", norm(F[:,1]))
    end

    x_old[:] = x
    f_old[:] = f
    return xnew

end