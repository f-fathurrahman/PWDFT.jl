using BenchmarkTools
using CuArrays

using PWDFT

function test_01()
    pw = PWGrid(15.0, gen_lattice_fcc(5.0))
    println(pw)

    ik = 1
    k = [0.0, 0.0, 0.0]

    Ngw = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]

    #Gk2 = zeros(ComplexF64,Ngw)
    psi = rand(ComplexF64,Ngw)
    
    #@views Gk2 = (pw.gvec.G[:,idx_gw2g] .+ k[:])
    Gk2 = transpose( sum( mapslices( x -> x.^2, pw.gvec.G[:,idx_gw2g] .+ k[:], dims=[2] ), dims=1 ) )

    println(size(Gk2))

    Kpsi = psi .* Gk2

    println(size(psi))
    #println(size(Kpsi))

end

test_01()