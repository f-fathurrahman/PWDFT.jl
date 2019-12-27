using CuArrays
using CuArrays.CUSOLVER
using LinearAlgebra
using Printf

CuArrays.allowscalar(false)

function test_01()

    N = 300

    A = rand(ComplexF64, N, N)
    A = 0.5*( A + A' )

    d_A = CuArray(A)
    d_A = 0.5*( d_A + d_A' )
    d_A_orig = copy(d_A)

    @time eigen_result = eigen(A)
    @time d_W, d_V = CUSOLVER.heevd!('V', 'U', d_A)
    
    d_A = copy(d_A_orig)  # set d_A to its original values

    @time eigen_result = eigen(A)
    @time d_W, d_V = CUSOLVER.heevd!('V', 'U', d_A)

    W = collect(d_W)

    if N <= 20
        for i in 1:N
            @printf("%5d %18.10f %18.10f\n", i, eigen_result.values[i], W[i])
        end
    end

    println("Pass here")
end
test_01()