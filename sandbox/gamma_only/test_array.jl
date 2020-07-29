using Printf
using LinearAlgebra
using Random

include("ortho_gram_schmidt.jl")
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


function overlap_gamma( psi1::Array{ComplexF64,2}, psi2::Array{ComplexF64,2} )
    C = psi2' * psi1
    Nstates = size(psi1,2)
    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi2[1,ist] # psi is the array that is conj transposed in the orig expression
        v2g[ist] = psi1[1,ist]  # v2g is the array that will be conj transposed
    end
    return C + conj(C) - v1g*v2g'
end

function test_ortho_sqrt()

    Nstates = 4

    psi1 = zeros(ComplexF64,5,Nstates)
    psi1[:,1] = [1.0, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]
    psi1[:,2] = [2.0, 4.1 + im*3, 5.2 + im*7, 4.1 - im*3, 5.2 - im*7]
    psi1[:,3] = [3.0, 9.0 + im*2, 8.1 + im*3, 9.0 - im*2, 8.1 - im*3]
    psi1[:,4] = [4.0, 1.1 + im*3, 1.2 + im*7, 1.1 - im*3, 1.2 - im*7]

    psi1g = zeros(ComplexF64,3,Nstates)
    psi1g[:,1] = [1.0, 3.0 + im*2, 4.1 + im*3]
    psi1g[:,2] = [2.0, 4.1 + im*3, 5.2 + im*7]
    psi1g[:,3] = [3.0, 9.0 + im*2, 8.1 + im*3]
    psi1g[:,4] = [4.0, 1.1 + im*3, 1.2 + im*7]

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

    println("dot psi1       = ", dot(psi1,psi1))
    println("dot psi2g      = ", 2*dot(psi1g,psi1g))
    ss = 0.0 + 0.0*im
    for ist in 1:Nstates
        ss = ss + 2*dot(psi1g[:,ist],psi1g[:,ist]) #- conj(psi1g[1,ist])*psi1g[1,ist]
    end
    println("gamma only dot = ", ss)
    println("dot 1, 2:", 2*dot( psi1g[:,1],psi1g[:,2] ) - conj(psi1g[1,1])*psi1g[1,2] )
    println("dot 1, 3:", 2*dot( psi1g[:,1],psi1g[:,3] ) - conj(psi1g[1,1])*psi1g[1,3] )

    println( dot(psi1[2:3,1],psi1[2:3,1]) )
    println( dot(psi1[2:3,1],psi1[2:3,2]) )

    println( dot(psi1[4:5,1],psi1[4:5,1]) )
    println( dot(psi1[4:5,1],psi1[4:5,2]) )
end
#test_ortho_sqrt()

function ortho_check_gamma( psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        #c = 2*dot( psi[:,ist], psi[:,ist] ) - conj(psi[1,ist])*psi[1,ist]
        c = dot( psi[:,ist], psi[:,ist] )
        c = c + conj(c) - conj(psi[1,ist])*psi[1,ist]
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        #c = 2*dot( psi[:,ist], psi[:,1] ) #- conj(psi[1,ist])*psi[1,1]
        c = dot( psi[:,ist], psi[:,1] )
        c = c + conj(c) - conj(psi[1,ist])*psi[1,1]
        #println(conj(psi[1,ist])*psi[1,1])
        #println(psi[1,ist]*conj(psi[1,1]))
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end

function gen_two_psi(Nbasis, Nstates)

    psig = randn(ComplexF64,Nbasis,Nstates)
    for ist in 1:Nstates
        #psig[1,ist] = 0.0 + im*0.0
        psig[1,ist] = psig[1,ist] + conj(psig[1,ist])
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
    Nbasis = 5 # for the gamma only

    psi1, psi1g = gen_two_psi(Nbasis,Nstates)
    
    println("\nBefore orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", 2*dot(psi1g,psi1g))

    G = psi1' * psi1
    Udagger = inv(sqrt(G))
    psi1o = psi1*Udagger

    G = psi1g' * psi1g
    Udagger = inv(sqrt(G + conj(G)))
    psi1go = psi1g*Udagger

    println("\nAfter orthonormalize ortho_sqrt")

    display(psi1o); println()
    display(psi1go); println()

    println("dot psi1       = ", dot(psi1o,psi1o))
    println("dot psi2g      = ", 2*dot(psi1go,psi1go))
    # dot product between different states are defined quite differently to have
    # result of zeros.
    for ist in 1:Nstates
        c1 = dot( psi1go[:,2], psi1go[:,ist] )
        c2 = dot( psi1go[:,ist], psi1go[:,2] )
        println("dot ist:", c1 + c2 )
    end

end
#test_ortho_sqrt2()

function test_ortho_GS_gamma2()
    
    Nstates = 4
    Nbasis = 6 # for the gamma only

    psi1, psi1g = gen_two_psi(Nbasis,Nstates)
    
    println("\nBefore orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", 2*dot(psi1g,psi1g))

    #psi1go = ortho_GS_gamma(psi1g)
    ortho_GS_gamma!(psi1g)

    println("\nAfter orthonormalize ortho_GS_gamma")

    #display(psi1o); println()
    display(psi1g); println()

    #println("dot psi1       = ", dot(psi1o,psi1o))
    println("dot psi2g      = ", 2*dot(psi1g,psi1g))
    for ist in 1:Nstates
        println("dot ist:", 2*dot( psi1g[:,ist],psi1g[:,4] ) )
    end

end
#test_ortho_GS_gamma2()

function dot_BlochWavefuncGamma( v1::Array{ComplexF64,2}, v2::Array{ComplexF64,2} )
    
    c = dot(v1, v2)
    s = c + conj(c)
    
    Nstates = size(v1,2)
    ds = 0.0 + im*0.0
    for ist in 1:Nstates
        ds = ds + conj(v1[1,ist]) * v2[1,ist]
    end
    return s - ds

end

function test_matrix2()
    
    Random.seed!(1234)

    Nstates = 3
    Nbasis = 4 # for the gamma only

    psi1, psi1g = gen_two_psi(Nbasis,Nstates)
    psi2, psi2g = gen_two_psi(Nbasis,Nstates)

    println("\nBefore orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot_BlochWavefuncGamma(psi1g,psi1g))

    #println("dot psi2  = ", dot(psi2,psi2))
    #println("dot psi2g = ", dot_BlochWavefuncGamma(psi2g,psi2g))

    ortho_gram_schmidt!(psi1)
    ortho_GS_gamma!(psi1g)

    #ortho_gram_schmidt!(psi2)
    #ortho_GS_gamma!(psi2g)

    println("\nAfter orthonormalize ortho_GS_gamma")

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot_BlochWavefuncGamma(psi1g,psi1g))

    println("\npsi1")
    display(psi1); println()

    println("\npsi1g")
    display(psi1g); println()


    U1 = psi1' * psi1

    C = psi1g' * psi1g
    U1g = C + conj(C)

    v1g = zeros(Float64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = real(psi1g[1,ist])
    end

    println()
    println("v1g = ")
    println(v1g)

    println()
    println("U1 and U1g")
    display(U1); println()
    display(U1g); println()

    uu = v1g*v1g'
    display(uu); println()

    display(U1g - uu); println()

end
#test_matrix2()


function test_matrix3()
    
    Random.seed!(1234)

    Nstates = 3
    Nbasis = 4 # for the gamma only

    psi1, psi1g = gen_two_psi(Nbasis,Nstates)
    psi2, psi2g = gen_two_psi(Nbasis,Nstates)

    println("\nBefore orthonormalize")

    display(psi1); println()
    display(psi1g); println()

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot_BlochWavefuncGamma(psi1g,psi1g))

    println("dot psi2  = ", dot(psi2,psi2))
    println("dot psi2g = ", dot_BlochWavefuncGamma(psi2g,psi2g))

    ortho_gram_schmidt!(psi1)
    ortho_GS_gamma!(psi1g)

    ortho_gram_schmidt!(psi2)
    ortho_GS_gamma!(psi2g)

    println("\nAfter orthonormalize ortho_GS_gamma")

    println("dot psi1  = ", dot(psi1,psi1))
    println("dot psi1g = ", dot_BlochWavefuncGamma(psi1g,psi1g))

    ortho_check_gamma(psi1g)

    println("dot psi1, psi2   = ", dot(psi1,psi2))
    println("dot psi1g, psi2g = ", dot_BlochWavefuncGamma(psi1g,psi2g))

    println("\npsi1")
    display(psi1); println()

    println("\npsi1g")
    display(psi1g); println()

    println("\npsi2")
    display(psi2); println()

    println("\npsi2g")
    display(psi2g); println()

    U12 = psi1' * psi2

    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi1g[1,ist]
        v2g[ist] = psi2g[1,ist]
    end

    println()
    println("U12 and U12")
    display(U12); println()
    #display(U12g); println()

    uu = v1g*v2g'
    #display(uu); println()

    U12g = U12g - uu
    display(U12g); println()

    println()
    println("Test multiply with U")

    println("\npsi1*U12")
    display(psi1*U12); println()

    println("\npsi1g")
    display(psi1g*U12g); println()

end
test_matrix3()

function constrain_search_dir( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    dnew = d - psi * ( psi' * d )
    return dnew
end

function constrain_search_dir!( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    d[:] = d - psi * ( psi' * d )
    return
end

function constrain_search_dir_gamma( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    C = psi' * d
    Nstates = size(d,2)
    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi[1,ist] # psi is the array that is conj transposed in the orig expression
        v2g[ist] = d[1,ist]  # v2g is the array that will be conj transposed
    end
    return d - psi * ( C + conj(C) - v1g*v2g' )
end

function constrain_search_dir_gamma!( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    C = psi' * d
    Nstates = size(d,2)
    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi[1,ist] # psi is the array that is conj transposed in the orig expression
        v2g[ist] = d[1,ist]  # v2g is the array that will be conj transposed
    end
    d[:] = d - psi * ( C + conj(C) - v1g*v2g' )
    return
end