using Printf
using LinearAlgebra

include("ortho_GS_gamma.jl")

function test_vector()
    
    # the DC component (first element) is real
    
    v1 = [3.1, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]

    v1g = [3.1, 3.0 + im*2, 4.1 + im*3] # gamma only

    println("dot v1  = ", dot(v1,v1))
    println("dot v1g = ", dot(v1g,v1g))
    println("gamma only dot = ", 2*dot(v1g,v1g) - conj(v1g[1])*v1g[1])
end
#test_vector()

function test_ortho_GS_gamma()

    Nstates = 4

    psi1 = zeros(ComplexF64,5,Nstates)
    psi1[:,1] = [0.0, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]
    psi1[:,2] = [0.0, 4.1 + im*3, 5.2 + im*7, 4.1 - im*3, 5.2 - im*7]
    psi1[:,3] = [0.0, 9.0 + im*2, 8.1 + im*3, 9.0 - im*2, 8.1 - im*3]
    psi1[:,4] = [0.0, 1.1 + im*3, 1.2 + im*7, 1.1 - im*3, 1.2 - im*7]

    psi1g = zeros(ComplexF64,3,Nstates)
    psi1g[:,1] = [0.0, 3.0 + im*2, 4.1 + im*3]
    psi1g[:,2] = [0.0, 4.1 + im*3, 5.2 + im*7]
    psi1g[:,3] = [0.0, 9.0 + im*2, 8.1 + im*3]
    psi1g[:,4] = [0.0, 1.1 + im*3, 1.2 + im*7]

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - psi1g[1,ist]*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)

    ortho_GS_gamma!(psi1g)

    println("dot psi1g = ", dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - conj(psi1g[1,ist])*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)

    display(psi1g); println()

end
#test_ortho_GS_gamma()

function test_ortho_sqrt()

    Nstates = 4

    psi1 = zeros(ComplexF64,5,Nstates)
    psi1[:,1] = [0.0, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]
    psi1[:,2] = [0.0, 4.1 + im*3, 5.2 + im*7, 4.1 - im*3, 5.2 - im*7]
    psi1[:,3] = [0.0, 9.0 + im*2, 8.1 + im*3, 9.0 - im*2, 8.1 - im*3]
    psi1[:,4] = [0.0, 1.1 + im*3, 1.2 + im*7, 1.1 - im*3, 1.2 - im*7]

    psi1g = zeros(ComplexF64,3,Nstates)
    psi1g[:,1] = [0.0, 3.0 + im*2, 4.1 + im*3]
    psi1g[:,2] = [0.0, 4.1 + im*3, 5.2 + im*7]
    psi1g[:,3] = [0.0, 9.0 + im*2, 8.1 + im*3]
    psi1g[:,4] = [0.0, 1.1 + im*3, 1.2 + im*7]

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - psi1g[1,ist]*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)

    G = psi1' * psi1
    println("G  = ")
    display(G); println()
    Udagger = inv(sqrt(G))
    psi1 = psi1*Udagger

    G = psi1g' * psi1g
    println("G  = ")
    display(G); println()
    println("conj(G) = ")
    display(conj.(G)); println()
    #display(G + G'); println()
    Udagger = inv(sqrt(G + conj(G)))
    psi1g = psi1g*Udagger

    println("\nAfter orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi2g = ", 2*dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) #- conj(psi1g[1,ist])*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)
    println("dot 1, 2:", 2*dot( psi1[:,1],psi1[:,2] ) )
    println("dot 1, 3:", 2*dot( psi1[:,1],psi1[:,3] ) )
end
#test_ortho_sqrt()

function gen_two_psi(Nbasis, Nstates)

    psig = randn(ComplexF64,Nbasis,Nstates)
    for ist in 1:Nstates
        psig[1,ist] = 0.0 + im*0.0
    end
    Nbasis2 = 2*Nbasis - 1
    psi = zeros(ComplexF64,Nbasis2,Nstates)
    for ist = 1:Nstates
        for ig in 1:Nbasis
            psi[ig,ist] = psig[ig,ist]
        end
        for ig in 2:Nbasis
            psi[Nbasis+ig-1,ist] = conj(psig[ig,ist])
        end
    end
    return psi, psig
end

function test_ortho_sqrt2()
    
    Nstates = 4
    Nbasis = 3 # for the gamma only

    psi1, psi1g = gen_two_psi(Nbasis,Nstates)
    
    println("\nBefore orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - psi1g[1,ist]*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)

    G = psi1' * psi1
    Udagger = inv(sqrt(G))
    psi1 = psi1*Udagger

    G = psi1g' * psi1g
    Udagger = inv(sqrt(G + conj(G)))
    psi1g = psi1g*Udagger

    println("\nAfter orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1       = ", dot(psi1,psi1))
    println("dot psi2g      = ", 2*dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) #- conj(psi1g[1,ist])*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)
    println("dot 1, 2:", 2*dot( psi1g[:,1],psi1g[:,2] ) )
    println("dot 1, 3:", 2*dot( psi1g[:,1],psi1g[:,3] ) )

end
test_ortho_sqrt2()