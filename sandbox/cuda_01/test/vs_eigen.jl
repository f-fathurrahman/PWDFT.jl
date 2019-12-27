using CuArrays
using CuArrays.CUSOLVER
using LinearAlgebra
using Printf

include("utils_CuArrays.jl")

function main( ; N=5)

    U = rand(ComplexF64,N,N)
    U = 0.5*(U + U')

    d_U = CuArray(U)

    println("N = ", N)

    @printf("Using CPU     : ")
    F = eigen(U)
    @time F = eigen(U)

    @printf("Using CUSOLVER: ")
    d_V, d_W = eigen(d_U)
    @time d_V, d_W = eigen(d_U)

    if N <= 5
        println(F.values)
        println(d_V)
    end
end

main(N=10)
main(N=50)
main(N=100)
main(N=200)
main(N=500)
main(N=1000)
