using Printf
using LinearAlgebra
using PWDFT

include("smearing.jl")

function test_01()
    evals = [1.0, 2.0, 3.0, 4.0]
    E_f = 2.5
    kT = 0.1
    
    println( smear_fermi.( evals, E_f, kT ) )

    x = (E_f .- evals)/kT
    println( PWDFT.wgauss.(x) )
end

#test_01()

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

test_03()