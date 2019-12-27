include("PWDFT_cuda.jl")
include("cu_Ylm_real.jl")

function kernel_call_Ylm_real( l, m, R1, R2, R3, vout )

    idx = ( blockIdx().x - 1 )*blockIdx().x + threadIdx().x

    N = length(R1)

    if idx <= N
        vout[idx] = cu_Ylm_real( l, m, R1[idx], R2[idx], R3[idx] )
    end

    return
end

function test_call_Ylm_real()
    Ndata = 100_000

    Nthreads = 256
    Nblocks = ceil(Int64, Ndata/Nthreads)

    R1 = CuArrays.rand(Float64, Ndata)
    R2 = CuArrays.rand(Float64, Ndata)
    R3 = CuArrays.rand(Float64, Ndata)
    vout = CuArrays.zeros(Float64, Ndata)

    for l in 0:3
        for m in -l:l
            @printf("l = %d m = %d\n", l, m)
            @cuda threads=Nthreads blocks=Nblocks kernel_call_Ylm_real( l, m, R1, R2, R3, vout )
        end
        @printf("\n")
    end

    println("Pass here")
end

test_call_Ylm_real()
