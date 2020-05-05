using Printf
using LinearAlgebra

include("ortho_GS_gamma_only.jl")

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

    ortho_GS_gamma_only!(psi1g)

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
end
test_ortho_sqrt()
