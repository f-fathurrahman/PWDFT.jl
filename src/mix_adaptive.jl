struct AdaptiveLinearMixer
    beta0::Float64
    betamax::Float64
    beta::Vector{Float64}
    f::Vector{Float64}
end

function AdaptiveLinearMixer(Rhoe::Array{Float64}, beta0; betamax=0.8)
    return AdaptiveLinearMixer(
        beta0,
        betamax,
        beta0*ones(Float64, length(Rhoe)),
        zeros(Float64, length(Rhoe))
    )
end

function do_mix!(
    mixer::AdaptiveLinearMixer,
    deltain, deltaout_,
    iterSCF::Int64
)
    # mixed quantity is deltain
    mix_adaptive!( deltain, deltaout_,
        mixer.beta0, mixer.beta, mixer.f, betamax=mixer.betamax
    )
    return
end


function mix_adaptive!( mu, nu, beta0::Float64, beta, f; betamax=0.8 )

    Npts = length(mu)
    
    for i in 1:Npts
        
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