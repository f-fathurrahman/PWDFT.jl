mutable struct BlochWavefuncGamma
    data::Vector{Array{ComplexF64,2}}
end

function BlochWavefuncGamma(psi::Array{ComplexF64,2})
    return BlochWavefuncGamma([psi])
end

# This assumes that the DC components are set to zeros
import LinearAlgebra: dot
function dot( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    
    #return 2*dot(v1.data, v2.data)
    
    c = dot(v1.data, v2.data)
    return c + conj(c)

    # alternative:
    # c = dot(v1.data, v2.data)
    # res = c + conj(c)
    # should work for ortho_sqrt and ortho_GS_gamma
end

import Base: length
function length(v::BlochWavefuncGamma)
    return length(v.data)
end

#
# Operators
#

import Base: +
function +( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    return BlochWavefuncGamma(v1.data + v2.data)
end

import Base: -
function -( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    return BlochWavefuncGamma(v1.data - v2.data)
end



function zeros_BlochWavefuncGamma( Ham::HamiltonianGamma )
    return zeros_BlochWavefuncGamma(Ham.pw.gvecw.Ngw, Ham.electrons.Nstates, Nspin=Ham.electrons.Nspin)
end


function zeros_BlochWavefuncGamma( Nbasis::Int64, Nstates::Int64; Nspin=1 )
    data = Vector{Array{ComplexF64,2}}(undef,Nspin)
    for ispin in 1:Nspin
        data[ispin] = zeros(ComplexF64,Nbasis,Nstates)
    end
    return BlochWavefuncGamma(data)
end




# random BlochWavefuncGamma from randn

function randn_BlochWavefuncGamma( Ham::HamiltonianGamma )
    return randn_BlochWavefuncGamma(Ham.pw.gvecw.Ngw, Ham.electrons.Nstates, Nspin=Ham.electrons.Nspin)
end

function randn_BlochWavefuncGamma( Nbasis::Int64, Nstates::Int64; Nspin=1 )

    data = Vector{Array{ComplexF64,2}}(undef,Nspin)
    
    G = zeros(ComplexF64,Nstates,Nstates)
    Udagger = zeros(ComplexF64,Nstates,Nstates)

    for ispin in 1:Nspin
        data[ispin] = randn(ComplexF64,Nbasis,Nstates)
        for ist in 1:Nstates
            data[ispin][1,ist] = 0.0 + im*0.0 # Don't forget to set the DC component to zero
        end
        G[:] = data[ispin]' * data[ispin]
        Udagger[:] = inv( sqrt(G + conj(G)) )
        data[ispin] = data[ispin]*Udagger
        #ortho_GS_gamma!(data[ispin])
    end

    return BlochWavefuncGamma(data)
end

function ortho_gram_schmidt!( psis::BlochWavefuncGamma )
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_GS_gamma!( psis.data[ispin] )
    end
    return
end

function ortho_check( psis::BlochWavefuncGamma )
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_check_gamma( psis.data[ispin] )
    end
    return
end

function ortho_check_gamma( psis::BlochWavefuncGamma )
    Nspin = length(psis)
    for ispin in 1:Nspin
        ortho_check_gamma( psis.data[ispin] )
    end
    return
end

function ortho_check_gamma( psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], psi[:,ist] )
        c = c + conj(c)
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        c = dot( psi[:,ist], psi[:,1] )
        c = c + conj(c)
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end


#
# Old definition, the DC component is not zero.
#
#function dot(v1::BlochWavefuncGamma, v2::BlochWavefuncGamma)
#    res = 2*dot(v1.data, v2.data)
#    for ist in 1:size(v1.data,2)
#        res = res - conj(v1.data[1,ist]) * v2.data[1,ist]
#    end
#    return res
#end