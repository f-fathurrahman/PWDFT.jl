using Printf
using LinearAlgebra
using PWDFT

include("smearing.jl")

function test_01()
    evals = [1.0, 2.0, 3.0, 4.0]
    Nstates = size(evals,1)
    E_f = 2.1
    kT = 0.1
    x = (E_f .- evals)/kT # for wgauss and w1gauss

    println( "smear        = ", smear_fermi.( evals, E_f, kT ) )
    println( "smear (v2)   = ", PWDFT.wgauss.(x) )

    x1 = -smear_fermi_entropy.( evals, E_f, kT ) 
    x2 = PWDFT.w1gauss.(x)
    println( "entropy      = ", x1)
    println( "entropy (v2) = ", x2)
    println( "ratio        = ", x2./x1)


    wks = [2.0]
    #mTS = calc_entropy( [1.0], kT, reshape(evals,(4,1)), E_f, 1 )
    mTS = 0.0
    println("\nOriginal")
    for ist = 1:Nstates
        xx = wks[1]*kT*PWDFT.w1gauss( (E_f - evals[ist])/kT )
        mTS = mTS + xx
        @printf("xx = %18.10f\n", xx)
    end
    println("mTS = ", mTS)

    mTS = 0.0
    println("\nNew")
    for ist = 1:Nstates
        xx = -wks[1]*(2*kT)*smear_fermi_entropy( evals[ist], E_f, kT )
        mTS = mTS + xx
        @printf("xx = %18.10f\n", xx)
    end
    println("mTS = ", mTS)
end
test_01()

function test_02()
    evals = [1.0, 2.0, 3.0, 4.0]
    E_f = 2.1
    kT = 0.1

    Nstates = length(evals)
    v = rand(ComplexF64,Nstates)
    ∇F = diagm( 0 => v )
    display(∇F); println()

    ∇ϵ = grad_smear( smear_fermi, smear_fermi_prime, evals, E_f, kT, ∇F )
    display(∇ϵ); println()
end

#test_02()

function test_03()
    evals = [1.0, 2.0, 3.0, 4.0]
    E_f = 2.1
    kT = 0.1

    Nstates = length(evals)
    fprime = zeros(Nstates)
    Focc = zeros(Nstates)
    for ist in 1:Nstates
        Focc[ist] = smear_fermi( evals[ist], E_f, kT )
        fprime[ist] = smear_fermi_prime( evals[ist], E_f, kT )
        df = Focc[ist]*(1.0 - Focc[ist])/kT
        @printf("%3d %18.10f %18.10f %18.10f\n", ist, Focc[ist], fprime[ist], df)
    end
end
#test_03()