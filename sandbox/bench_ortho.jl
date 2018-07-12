using LinearAlgebra
using PWDFT
using BenchmarkTools

const Ngw = 2000
const Nstates = 20

function using_sqrt()
    psi = ortho_sqrt(rand(ComplexF64,Ngw,Nstates))
end

function using_gram_schmidt()
    psi = ortho_gram_schmidt(rand(ComplexF64,Ngw,Nstates))
end

@btime using_sqrt()
@btime using_gram_schmidt()

#=

Size: (1000,10)
  509.231 μs (30 allocations: 344.59 KiB): using_sqrt
  1.886 ms (479 allocations: 5.08 MiB): using_gram_schmidt

Size: (2000,20)
  2.535 ms (30 allocations: 1.30 MiB)
  14.198 ms (3074 allocations: 38.56 MiB)
=#