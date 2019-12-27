using Test
using Random
using Printf

using PWDFT
using PWDFT_cuda

function main()

    pw = CuPWGrid(15.0, gen_lattice_fcc(10.0))

    psiks = zeros_CuBlochWavefunc(pw, 10, 1)

    println(typeof(psiks))

    psiks2 = rand_CuBlochWavefunc(pw, 10, 1)
    ortho_check(psiks2[1])

    println("Pass here")

end

main()