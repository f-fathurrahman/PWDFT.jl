
# f_n = \left( 1 + \exp\left[ \epsilon_{n} - \mu \right]/kT \right)^{-1}
# kT = smearing width
function smear_FD( ev::Array{Float64,1}, efermi::Float64, kT::Float64; Nspin=1 )
    #
    x = (ev[:] - efermi)/kT
    #
    Nstates = length(ev)
    f = zeros(Nstates)
    if Nspin == 2
        for ist = 1:Nstates
            f[ist] = 1.0/( 1.0 + exp(x[ist]))
        end
    else
        for ist = 1:Nstates
            f[ist] = 2.0/( 1.0 + exp(x[ist]))
        end        
    end
    return f
end


function smear_FD( ev::Float64, efermi::Float64, kT::Float64; Nspin=1 )
    x = (ev - efermi)/kT
    if Nspin == 2
        f = 1.0/( 1.0 + exp(x))
    else
        f = 2.0/( 1.0 + exp(x))
    end
    return f
end
