using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

include("PWGridGamma.jl")

function main()
    
    Random.seed!(1234)

    LatVecs = gen_lattice_sc(6.0)

    pw = PWGrid(5.0, LatVecs)

    pw_gamma = PWGridGamma(5.0, LatVecs)

    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    G = pw_gamma.gvec.G
    G2 = pw_gamma.gvec.G2
    for ig = 1:4
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end

    # Ordinary G
    println()
    for ig in 1:4
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig,
            pw.gvec.G[1,ig], pw.gvec.G[2,ig], pw.gvec.G[3,ig], pw.gvec.G2[ig])
    end

end

main()