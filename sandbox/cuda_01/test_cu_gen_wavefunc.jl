using Printf
using FFTW
using LinearAlgebra

using CuArrays

using PWDFT

include("CuPWGrid.jl")

include("cu_types_aliases.jl")

include("cu_ortho_gram_schmidt.jl")
include("cu_gen_wavefunc.jl")
include("cu_ortho_check.jl")

function main()

    pw = CuPWGrid(15.0, gen_lattice_fcc(10.0))

    psiks = zeros_CuBlochWavefunc(pw, 10, 3)

    println(typeof(psiks))

    psiks2 = rand_CuBlochWavefunc(pw, 10, 3)
    ortho_check(psiks2[1])

    println("Pass here")

end

main()