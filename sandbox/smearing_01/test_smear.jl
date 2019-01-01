using Printf
using PWDFT

#
# Smear FD
#
function smear_FD_v2( ev::Array{Float64,1}, efermi::Float64, kT::Float64; Nspin=1 )
    #
    x = (ev[:] .- efermi)/kT
    #
    Nstates = length(ev)
    f = zeros(Nstates)
    if Nspin==1
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

#=
This function need to be be generalized to accomodate Nkspin and ikspin
=#
function test_main(kT::Float64)
    Nstates = 5
    Nspin = 2
    
    Focc    = zeros(Nstates,Nspin)
    Focc_v2 = zeros(Nstates,Nspin)
    
    ev = Array{Float64}(undef,Nstates,Nspin)
    ev[:,1] = [-2.4, -1.0, -0.5, -0.2, -0.19]
    ev[:,2] = [-2.4, -1.0, -0.5, -0.2, -0.19]
    
    wk = ones(Nstates,1)
    
    efermi = 0.5*(ev[3] + ev[4])
    
    println("\nkT = ", kT)
    
    println("\nNspin = ", Nspin)
    Focc[:,1] = smear_FD(ev[:,1], efermi, kT, Nspin=Nspin)
    Focc_v2[:,2] = smear_FD_v2(ev[:,2], efermi, kT, Nspin=Nspin)
    for ist = 1:Nstates
        @printf("%18.10f %18.10f %18.10f\n", ev[ist,1], Focc[ist,1], Focc_v2[ist,1])
    end
    @printf("sum(Focc) = %18.10f\n", sum(Focc))
    @printf("Entropy (-TS) = %18.10f\n", calc_entropy(Focc, wk[:,1], kT, Nspin=Nspin))
end


test_main(0.001)
test_main(0.01)
test_main(0.1)
