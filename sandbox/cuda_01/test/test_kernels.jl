using Test
using CUDAnative
using CuArrays
using SpecialFunctions

function kernel_multiply( vin, vout )
    idx = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    N = length(vin)
    if idx <= N
        @inbounds vout[idx] = 2.0*vin[idx]
    end
    return
end

function test_multiply()
    Ndata = 100_000

    Nthreads = 256
    Nblocks = ceil(Int64, Ndata/Nthreads)

    vin = CuArrays.rand(Float64, Ndata)
    vout = CuArrays.zeros(Float64, Ndata)

    @cuda threads=Nthreads blocks=Nblocks kernel_multiply( vin, vout )

    println("Pass here")
end

#test_multiply()


function kernel_call_myfunc1( C1, C2, C3, vin, vout )
    idx = ( blockIdx().x - 1 )*blockIdx().x + threadIdx().x
    N = length(vin)
    if idx <= N
        @inbounds vout[idx] = myfunc1( C1, C2, C3, vin[idx] )
    end
    return
end

function myfunc1( C1, C2, C3, x )
    return C1*x*x + C2*x + C3
end

function test_call_myfunc1()
    Ndata = 100_000

    Nthreads = 256
    Nblocks = ceil(Int64, Ndata/Nthreads)

    vin = CuArrays.rand(Float64, Ndata)
    vout = CuArrays.zeros(Float64, Ndata)

    C1 = 1.0
    C2 = 2.0
    C3 = 3.0
    @cuda threads=Nthreads blocks=Nblocks kernel_call_myfunc1( C1, C2, C3, vin, vout )

    vin_ref = collect(vin)
    vout_ref = myfunc1.(C1, C2, C3, vin_ref)

    vout_res = collect(vout)

    @test vout_ref ≈ vout_ref

    println("Pass here call myfunc1")
end

test_call_myfunc1()
