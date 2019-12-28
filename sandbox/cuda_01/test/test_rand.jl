using Random
using CuArrays

Random.seed!(1234)
CuArrays.CURAND.seed!(1234)

function gen_rand()
    v = rand(ComplexF64,2,2)
    return v
end

function cu_gen_rand()
    v = CuArrays.rand(ComplexF64,2,2)
    return v
end

function main()
    v_cpu = gen_rand()
    v_gpu = gen_rand()

    println(v_cpu)
    println(v_gpu)
end

main()