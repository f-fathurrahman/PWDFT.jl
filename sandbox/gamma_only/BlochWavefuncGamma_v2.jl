# Using separate array for dc components

mutable struct BlochWavefuncGamma
    dc::Vector{Vector{Float64}}  # FIXME use real instead of complex?
    data::Vector{Array{ComplexF64,2}}
end

# Constructor, from an array to BlochWavefuncGamma (Nspin=1)
function BlochWavefuncGamma(psi_::Array{ComplexF64,2})
    psi = copy(psi_) # do not modify argument
    Nstates = size(psi,2)
    dc = zeros(Float64,Nstates)
    for ist in 1:Nstates
        dc[ist] = real(psi[1,ist])
        psi[1,ist] = 0.0 + im*0.0
    end
    return BlochWavefuncGamma([dc], [psi])
end


# version 2
function dot_BlochWavefuncGamma( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    res = 0.0
    Nspin = length(v1.data)
    for ispin in 1:Nspin
        c = dot(v1.data[ispin], v2.data[ispin])
        res = res + c + conj(c) + dot(v1.dc[ispin], v2.dc[ispin])
    end
    return res
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
    return BlochWavefuncGamma(v1.dc + v2.dc, v1.data + v2.data)
end

import Base: -
function -( v1::BlochWavefuncGamma, v2::BlochWavefuncGamma )
    return BlochWavefuncGamma(v1.dc - v2.dc, v1.data - v2.data)
end



function zeros_BlochWavefuncGamma( Ham::HamiltonianGamma )
    return zeros_BlochWavefuncGamma(Ham.pw.gvecw.Ngw, Ham.electrons.Nstates, Nspin=Ham.electrons.Nspin)
end

function zeros_BlochWavefuncGamma( Nbasis::Int64, Nstates::Int64; Nspin=1 )
    dc = Vector{Vector{Float64}}(undef,Nspin)
    data = Vector{Array{ComplexF64,2}}(undef,Nspin)
    for i in 1:Nspin
        dc[i] = zeros(Float64,Nstates)
        data[i] = zeros(ComplexF64,Nbasis,Nstates)
    end
    return BlochWavefuncGamma(dc, data)
end


# random BlochWavefuncGamma from randn

function randn_BlochWavefuncGamma( Ham::HamiltonianGamma )
    return randn_BlochWavefuncGamma(Ham.pw.gvecw.Ngw, Ham.electrons.Nstates, Nspin=Ham.electrons.Nspin)
end

function randn_BlochWavefuncGamma( Nbasis::Int64, Nstates::Int64; Nspin=1 )
    dc = Vector{Vector{Float64}}(undef,Nspin)
    data = Vector{Array{ComplexF64,2}}(undef,Nspin)
    for ispin in 1:Nspin
        data[ispin] = randn(ComplexF64,Nbasis,Nstates)
        for ist in 1:Nstates
            dc[ispin][ist] = real(data[ispin][1,ist])
            data[ispin][1,ist] = 0.0 + im*0.0
        end
        ortho_sqrt_gamma!(dc[ispin], data[ispin])
    end
    return BlochWavefuncGamma(dc, data)
end


function overlap_gamma(
    dc1::Vector{Float64}, psi1::Array{ComplexF64,2},
    dc2::Vector{Float64}, psi2::Array{ComplexF64,2}
)
    C = psi2' * psi1
    return C + conj(C) + dc1*dc2'
end

function overlap_gamma( i::Int64, psis1::BlochWavefuncGamma, psis2::BlochWavefuncGamma )
    return overlap_gamma( psis1.dc[i], psis1.data[i],
                          psis2.dc[i], psis2.data[i] )
end

function ortho_sqrt_gamma!( psis::BlochWavefuncGamma )
    Nspin = size(psis.data)
    for i in 1:Nspin
        ortho_sqrt_gamma(psis.dc[i], psis.data[i])
    end
    return
end

function ortho_sqrt_gamma!( dc::Vector{Floa64}, psi::Array{ComplexF64,2} )
    Udagger = inv(sqrt(overlap_gamma(psi, psi)))
    psi[:,:] = psi[:,:]*Udagger
    dc[:] = Udagger*dc[:]
    return
end

# FIXME:
# remove _gamma postfix for ortho_sqrt, differentiated by number of arguments
# and arguments' type.


#
# Using raw array: dc components are the first element of each columns
#

#
function dot_BlochWavefuncGamma( v1::Array{ComplexF64,2}, v2::Array{ComplexF64,2} )
    c = dot(v1, v2)
    s = c + conj(c)
    Nstates = size(v1,2)
    ds = 0.0 + im*0.0
    for ist in 1:Nstates
        ds = ds + conj(v1[1,ist]) * v2[1,ist]
    end
    return s - ds
end

# For calculating psi2' * psi1
# order of argument is importatn
function overlap_gamma( psi1::Array{ComplexF64,2}, psi2::Array{ComplexF64,2} )
    C = psi2' * psi1
    Nstates = size(psi1,2)
    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi2[1,ist]  # v1g is the array that is conj transposed in the orig expression
        v2g[ist] = psi1[1,ist]  # v2g is the array that will be conj transposed
    end
    return C + conj(C) - v1g*v2g'
end

function ortho_sqrt_gamma!( psi::Array{ComplexF64,2} )
    C = overlap_gamma(psi, psi)
    psi[:,:] = psi[:,:]*inv(sqrt(C))
    return
end

# Return new array
function ortho_sqrt_gamma( psi::Array{ComplexF64,2} )
    C = overlap_gamma(psi, psi)
    return psi[:,:]*inv(sqrt(C))
end