include("PWDFT_cuda.jl")
include("kernels_PsPot_GTH.jl")

function kernel_call_eval_proj_G( rrl, l, iproj, CellVolume, vin, vout )

    idx = ( blockIdx().x - 1 )*blockIdx().x + threadIdx().x

    N = length(vin)

    if idx <= N
        vout[idx] = cu_eval_proj_G( rrl, l, iproj, vin[idx], CellVolume )
    end

    return
end

function test_call_eval_proj_G()
    Ndata = 100_000

    Nthreads = 256
    Nblocks = ceil(Int64, Ndata/Nthreads)

    vin = abs.( CuArrays.rand(Float64, Ndata) )
    vout = CuArrays.zeros(Float64, Ndata)

    rrl = 1.0
    l = 2
    iproj = 1
    CellVolume = 3.0
    @cuda threads=Nthreads blocks=Nblocks kernel_call_eval_proj_G( rrl, l, iproj, CellVolume, vin, vout )

    println("Pass here")
end

test_call_eval_proj_G()
