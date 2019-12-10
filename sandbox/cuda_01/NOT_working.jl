# This is not working because C is not isbits
function kernel_call_myfunc2( C, vin, vout )
    idx = ( blockIdx().x - 1 )*blockIdx().x + threadIdx().x
    N = length(vin)
    if idx <= N
        @inbounds vout[idx] = myfunc2( C, vin[idx] )
    end
    return
end

function myfunc2( C::Array{Float64,1}, x )
    return C[1]*x*x + C[2]*x + C[3]
end

function test_call_myfunc2()
    Ndata = 100_000

    Nthreads = 256
    Nblocks = ceil(Int64, Ndata/Nthreads)

    vin = CuArrays.rand(Float64, Ndata)
    vout = CuArrays.zeros(Float64, Ndata)

    C = [1.0, 2.0, 3.0]
    @cuda threads=Nthreads blocks=Nblocks kernel_call_myfunc2( C, vin, vout )

    vin_ref = collect(vin)
    vout_ref = zeros(Float64, Ndata)
    for i in 1:Ndata
        vout_ref[i] = myfunc2(C, vin_ref[i])
    end

    vout_res = collect(vout)

    @test vout_ref ≈ vout_ref

    println("Pass here call myfunc2")
end

test_call_myfunc2()
