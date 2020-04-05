# x = (ϵ - E_f)/(2.0*kT)

function smear_fermi( ϵ, E_f::Float64, kT::Float64 )
   x = (ϵ - E_f)/(2.0*kT)
   return 0.5*(1.0 - tanh(x))
end

function smear_gauss( ϵ, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return 0.5*erfc(x)
end

function smear_cold( ϵ, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return 0.5*erfc(x+sqrt(0.5)) + exp( -(x + sqrt(0.5))^2 )/sqrt(2*pi)
end

function smear_fermi_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    f = 0.5*( 1 - tanh(x) )
    S = 0.0
    if f > 1e-300
        S = S - f*log(f)
    elseif (1-f) > 1e-300
        S = S - (1-f)*log(1-f)
    end
    return S
end

function smear_gauss_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return exp(-x^2) / sqrt(pi)
end

function smear_cold_entropy( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return exp(-(x+sqrt(0.5)^2)) * (1 + x*sqrt(2))/sqrt(pi)
end

function smear_fermi_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -0.25/(kT*cosh(x)^2)
end


function smear_gauss_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -exp(-x*x) / (2.0*sqrt(pi)*kT)
end

function smear_cold_prime( ϵ::Float64, E_f::Float64, kT::Float64 )
    x = (ϵ - E_f)/(2.0*kT)
    return -exp( -(x+sqrt(0.5))^2 ) * (2 + x*sqrt(2)) / (2*sqrt(pi)*kT)
end
