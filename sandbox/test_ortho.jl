using Random
using PWDFT

function using_gram_schmidt()
    Random.seed!(2345)
    psi = rand(1000,4) + im*rand(1000,4)
    #
    println("\nBefore ortho_gram_schmidt:")
    ortho_check(psi)
    #
    psi = ortho_gram_schmidt(psi)
    println("\nAfter ortho_gram_schmidt:")
    ortho_check(psi)
end

function using_sqrt()
    Random.seed!(2345)
    psi = rand(1000,4) + im*rand(1000,4)
    #
    println("\nBefore ortho_sqrt:")
    ortho_check(psi)
    #
    psi = ortho_sqrt(psi)
    println("\nAfter ortho_sqrt:")
    ortho_check(psi)
end

function test_ortho_sqrt()
    psi = rand(ComplexF64,Ngw,Nstates)
    ortho_sqrt!(psi)
    ortho_check(psi)
end

using_gram_schmidt()
using_sqrt()

test_ortho_sqrt()
