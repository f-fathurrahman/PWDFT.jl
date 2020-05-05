mutable struct BlochWavefuncGamma
    data::Array{ComplexF64,2}
end

import LinearAlgebra: dot
function dot( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    return 2*dot(v1.data, v2.data)
end

function randn_BlochWavefuncGamma( Ham::Hamiltonian )
    return randn_BlochWavefuncGamma(Ham.pw.gvec.Ng, Ham.electrons.Nstates)
end

function randn_BlochWavefuncGamma( Nbasis::Int64, Nstates::Int64 )
    data = randn(ComplexF64,Nbasis,Nstates)
    for ist in 1:Nstates
        data[1,ist] = 0.0 + im*0.0 # Don't forget to set the DC component to zero
    end
    #G = data' * data
    #Udagger = inv( sqrt(G + conj.(G)) )
    #data = data*Udagger
    ortho_GS_gamma!(data)
    return BlochWavefuncGamma(data)
end

function ortho_check( psi::BlochWavefuncGamma )
    ortho_check_gamma( psi.data )
    return
end

function ortho_check_gamma( psi::BlochWavefuncGamma )
    ortho_check_gamma( psi.data )
    return
end

function ortho_check_gamma( psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = 2*dot( psi[:,ist], psi[:,ist] )
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        c = 2*dot( psi[:,ist], psi[:,1] )
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end


#
# Old definition, the DC component is not zero, but a pure real number
#
#function dot(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
#    res = 2*dot(v1.data, v2.data)
#    for ist in 1:size(v1.data,2)
#        res = res - v1.data[1,ist] * v2.data[1,ist]
#    end
#    return res
#end