function mix_rpulay!( x, gx, beta, X, F, iter, MIXDIM )

    f = gx - x
    println("norm f = ", norm(f))

# iter = 1, active_dim = 1
# iter = 2, active_dim = 2
# ...
# iter = MIXDIM, active_dim = MIXDIM
# iter = MIXDIM+1, active_dim = 1
# iter = MIXDIM+2, active_dim = 2
# ...

    active_dim = mod(iter,MIXDIM)
    if active_dim == 0
        active_dim = 4
    end

    println("\nin mix_rpulay:")
    println("active_dim = ", active_dim)

    if iter > 1
        # Restart, use the last history
        if active_dim == 1
            X[:,1] = X[:,MIXDIM]
            F[:,1] = F[:,MIXDIM]
            X[:,2:MIXDIM] .= 0.0
            F[:,2:MIXDIM] .= 0.0
        else
            for i = active_dim:-1:2
                X[:,i] = X[:,i-1]
                F[:,i] = F[:,i-1]
                println(i, "norm X[:,i] = ", norm(X[:,i]))
                println(i, "norm F[:,i] = ", norm(F[:,i]))
            end
            X[:,1] = x - X[:,2]
            F[:,1] = f - F[:,2]
            println("norm X[:,1] = ", norm(X[:,1]))
            println("norm F[:,1] = ", norm(F[:,1]))
        end
        Xk = X[:,1:active_dim]
        Fk = F[:,1:active_dim]
        #if active_dim == 1
        #    println("Pass here 34 in mix_rpulay")
        #    cc = Fk'*Fk
        #    println("cc = ", cc)
        #end
        addv = (Xk + beta*Fk)*inv(Fk'*Fk)*(Fk'*f)
        xnew = x + beta*f - addv
        println("norm addv = ", norm(addv))
        println("norm xnew = ", norm(xnew))
    else
        xnew = x + beta*f
        X[:,1] = x
        F[:,1] = f
        println("norm X[:,1] = ", norm(X[:,1]))
        println("norm F[:,1] = ", norm(F[:,1]))
    end

    return xnew

end