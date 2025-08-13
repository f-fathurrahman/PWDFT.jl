# x = (ϵ - E_f)/(2.0*kT)

function smear_fermi( ϵ, E_f::Float64, kT::Float64 )
   x = (ϵ - E_f)/(2.0*kT)
   return 0.5*(1.0 - tanh(x))
end

function smear_fermi_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    f = 0.5*( 1 - tanh(x) )
    S = 0.0
    if f > 1e-300
        S = S - f*log(f)
    end
    if (1-f) > 1e-300
        S = S - (1-f)*log(1-f)
    end
    return S
end

function smear_fermi_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -0.25/(kT*cosh(x)^2)
end


#
# Gaussian smearing
#

function smear_gauss( ϵ, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return 0.5*erfc(x)
end

function smear_gauss_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return exp(-x^2) / sqrt(pi)
end

function smear_gauss_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -exp(-x*x) / (2.0*sqrt(pi)*kT)
end

#
# Cold smearing
#
function smear_cold( ϵ, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return 0.5*erfc(x+sqrt(0.5)) + exp( -(x + sqrt(0.5))^2 )/sqrt(2*pi)
end

function smear_cold_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return exp( -(x+sqrt(0.5))^2 ) * (1 + x*sqrt(2))/sqrt(pi)
end

function smear_cold_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -exp( -(x+sqrt(0.5))^2 ) * (2 + x*sqrt(2)) / (2*sqrt(pi)*kT)
end


#
# Gradient
#

# `smear` and `smear_prime` must match
function grad_smear(
    smear_func, smear_func_prime,
    evals::Array{Float64,1}, E_f::Float64, kT::Float64, ∇F
)
    ∇ϵ = copy(∇F)
    #
    Nrows = size(∇ϵ,1)
    Ncols = size(∇ϵ,2)
    for j in 1:Ncols, i in 1:Nrows
        dϵ = evals[i] - evals[j]
        if abs(dϵ) < 1e-6
            ∇ϵ[i,j] = ∇ϵ[i,j] * smear_func_prime( evals[i], E_f, kT )
        else
            ∇ϵ[i,j] = ∇ϵ[i,j] * ( smear_func(evals[i], E_f, kT) - smear_func(evals[j], E_f, kT) ) / dϵ
        end
    end
    return ∇ϵ
end


# Other versions

# From QE
# x = E_f - evals
# only n=-99 is implemented
function wgauss( ϵ, E_f::Float64, kT::Float64 )
    #
    x = (E_f - ϵ)/kT
    #
    MAXARG = 200.0
    if x < -MAXARG
        return 0.0
    elseif x > MAXARG
        return 1.0
    else
        return 1.0/(1.0 + exp(-x))
    end

end

# For calculating entropy, the sign is different from smear_*_entropy
function w1gauss(ϵ, E_f::Float64, kT::Float64 )
    #
    x = (E_f - ϵ)/kT
    # n = -99 other cases are not yet implemented
    if abs(x) <= 36.0
        f = 1.0/(1.0 + exp(-x))
        onemf = 1.0 - f
        return f*log(f) + onemf*log(onemf)
    else
        return 0.0
    end
end