using LinearAlgebra
using BenchmarkTools

function test1()
    dw = zeros(1000,1000)
    v = rand(1000)
    h = rand(1000)
    for i = 1:10
      dw += h*v'
    end
    #BLAS.gemm!('N', 'T', 1.0, h, v, 1.0, dw)
end

function test2()
    dw = zeros(1000,1000)
    v = rand(1000)
    h = rand(1000)
    for i = 1:10
        #dw += h*v'
        BLAS.gemm!('N', 'T', 1.0, h, v, 1.0, dw)
    end
end

@btime test1()
@btime test2()