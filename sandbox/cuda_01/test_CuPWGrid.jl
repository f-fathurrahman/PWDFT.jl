using CuArrays
using LinearAlgebra
using FFTW

using PWDFT

include("CuPWGrid.jl")

function main()
    pw = CuPWGrid( 15.0, gen_lattice_fcc(10.0) )
    println("Pass here")
end

main()