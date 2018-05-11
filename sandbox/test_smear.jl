using PWDFT

#
# Smear FD
#
function smear_FD_v2( ev::Array{Float64,1}, efermi::Float64, kT::Float64; is_spinpol=false )
    #
    x = (ev[:] - efermi)/kT
    #
    Nstates = length(ev)
    f = zeros(Nstates)
    if is_spinpol
        for ist = 1:Nstates
            #f[ist] = 1.0/( 1.0 + exp(x[ist]))
            f[ist] = 0.25/( cosh(x[ist]/2.0)^2 )
        end
    else
        for ist = 1:Nstates
            #f[ist] = 2.0/( 1.0 + exp(x[ist]))
            f[ist] = 0.5/( cosh(x[ist]/2.0)^2 )
        end        
    end
    return f
end

function test_main(kT::Float64)
    Nstates = 5
    ev = Array{Float64}(Nstates)
    ev = [-2.4, -1.0, -0.5, -0.2, -0.19]
    efermi = 0.5*(ev[3] + ev[4])
    println("\nkT = ", kT)
    for spinpol in [true,false]
        println("\nspinpol = ", spinpol)
        Focc = smear_FD(ev, efermi, kT, is_spinpol=spinpol)
        Focc_v2 = smear_FD_v2(ev, efermi, kT, is_spinpol=spinpol)
        for ist = 1:Nstates
            @printf("%18.10f %18.10f %18.10f\n", ev[ist], Focc[ist], Focc_v2[ist])
        end
        #@printf("sum(Focc) = %18.10f\n", sum(Focc))
        #@printf("Entropy (-TS) = %18.10f\n", calc_entropy(Focc, kT, is_spinpol=spinpol))
    end
end


#test_main(0.001)
test_main(0.01)
test_main(0.1)
