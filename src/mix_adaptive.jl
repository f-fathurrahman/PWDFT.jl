function mix_adaptive!( mu, nu, beta0::Float64, beta, f; betamax=0.8 )

    Npts = length(mu)
    
    for i = 1:Npts
        
        t = nu[i] - mu[i]

        if t*f[i] >= 0.0
            beta[i] = beta[i] + beta0
            if beta[i] > betamax
                beta[i] = betamax
            end
        else
            beta[i] = 0.5*( beta[i] + beta0 )
        end
        
        f[i] = t
        
        mu[i] = beta[i]*nu[i] + ( 1.0 - beta[i] )*mu[i]

    end

    return

end