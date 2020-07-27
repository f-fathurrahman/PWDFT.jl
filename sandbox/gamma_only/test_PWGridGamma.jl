using Printf
using LinearAlgebra
using FFTW

using PWDFT

include("PWGridGamma.jl")

function main(LatVecs)

    # The usual PWGrid
    pw = PWGrid(15.0, LatVecs)
    
    # Using Gamma-point trick
    pw_gamma = PWGridGamma(15.0, LatVecs)

    # Size comparison
    println("summarysize(pw_gamma) = ", Base.summarysize(pw_gamma))
    println("summarysize(pw)       = ", Base.summarysize(pw))

    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    println("G-vectors of pw (Gamma-point trick)")
    G = pw_gamma.gvec.G
    G2 = pw_gamma.gvec.G2
    for ig = 1:4
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end

    # Ordinary G
    println()
    println("G-vectors of pw (usual)")
    for ig in 1:10
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig,
            pw.gvec.G[1,ig], pw.gvec.G[2,ig], pw.gvec.G[3,ig], pw.gvec.G2[ig])
    end

    println("Ng (Gamma-point trick) = ", pw_gamma.gvec.Ng)
    println("Ng (usual)             = ", pw.gvec.Ng)

    println("Ngw (Gamma-point trick) = ", pw_gamma.gvecw.Ngw)
    println("Ngw (usual) = ", pw.gvecw.Ngw)

end
#main(gen_lattice_sc(16.0))
#main(gen_lattice_fcc(16.0))
main(diagm([10.0, 11.0, 12.0]))
