using Printf
using LinearAlgebra

using PWDFT

include("init_gvec_gamma.jl")

function main()
    
    LatVecs = gen_lattice_sc(5.0)

    pw = PWGrid(5.0, LatVecs)
    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    init_gvec_gamma(Ns, RecVecs, ecutrho)
end

main()