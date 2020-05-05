using Printf
using LinearAlgebra

function test_vector()
    
    # the DC component (first element) is real
    
    v1 = [3.1, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]

    v1g = [3.1, 3.0 + im*2, 4.1 + im*3] # gamma only

    println("dot v1  = ", dot(v1,v1))
    println("dot v1g = ", dot(v1g,v1g))
    println("gamma only dot = ", 2*dot(v1g,v1g) - v1g[1]*v1g[1])
end
#test_vector()

function test_matrix()

    Nstates = 4

    psi1 = zeros(ComplexF64,5,Nstates)
    psi1[:,1] = [3.1, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]
    psi1[:,2] = [2.1, 4.1 + im*3, 5.2 + im*7, 4.1 - im*3, 5.2 - im*7]
    psi1[:,3] = [5.1, 9.0 + im*2, 8.1 + im*3, 9.0 - im*2, 8.1 - im*3]
    psi1[:,4] = [6.1, 1.1 + im*3, 1.2 + im*7, 1.1 - im*3, 1.2 - im*7]

    psi1g = zeros(ComplexF64,3,Nstates)
    psi1g[:,1] = [3.1, 3.0 + im*2, 4.1 + im*3]
    psi1g[:,2] = [2.1, 4.1 + im*3, 5.2 + im*7]
    psi1g[:,3] = [5.1, 9.0 + im*2, 8.1 + im*3]
    psi1g[:,4] = [6.1, 1.1 + im*3, 1.2 + im*7]

    println("dot v1  = ", dot(psi1,psi1))
    println("dot v1g = ", dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - psi1g[1,ist]*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)

    p = psi1[2:end,:]
    G = p' * p
    display(G); println()

    p = psi1g[2:end,:]
    Gg = p' * p
    display(Gg); println()

    #Udagger = inv( sqrt( psi1' * psi1 ) )
    #psi1 = psi1*Udagger

    #Udagger = inv( sqrt( psi1g' * psi1g ) )
    #psi1g = psi1g*Udagger

    #println("\nAfter orthonormalize")
    #println("dot v1  = ", dot(psi1,psi1))
    #ss = 0.0 + 0.0*im
    #for ist in 1:Nstates
    #    ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) - psi1g[1,ist]*psi1g[1,ist]
    #end
    #println("gamma only dot = ", ss)
end

test_matrix()
