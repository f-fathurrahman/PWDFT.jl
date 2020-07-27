using Printf
using LinearAlgebra
using FFTW

using PWDFT

include("init_gvec_gamma_debug2.jl")

function main()
    
    #LatVecs = diagm([10.0, 11.526478000000001, 10.596309])

    #LatVecs = diagm([10.0, 10.0, 10.0])
    
    LatVecs = diagm([2.0, 3.0, 5.0])

    # The usual PWGrid
    pw = PWGrid(5.0, LatVecs)
    #println(pw)

    println("pw.Ns = ", pw.Ns)
    init_gvec_gamma( pw.Ns, pw.RecVecs, pw.ecutrho )

end

main()