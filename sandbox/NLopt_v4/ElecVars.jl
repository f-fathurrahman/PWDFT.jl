mutable struct SubspaceRotations
    prev::Vector{Matrix{ComplexF64}}
    prevC::Vector{Matrix{ComplexF64}}
    prevCinv::Vector{Matrix{ComplexF64}}
end

function SubspaceRotations(Nkspin::Int64, Nstates::Int64)
    prev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    prevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    prevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        prev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        prevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        prevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end
    return SubspaceRotations(prev, prevC, prevCinv)
end

mutable struct ElecGradient
    psiks::BlochWavefunc
    Haux::Vector{Matrix{ComplexF64}}
end

function ElecGradient(Ham)
    psiks = zeros_BlochWavefunc(Ham)
    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = zeros(ComplexF64,Nstates,Nstates)
    end
    return ElecGradient(psiks, Haux)
end

import Base: -
function -(e1::ElecGradient, e2::ElecGradient)
    return ElecGradient( e1.psiks .- e2.psiks, e1.Haux .- e2.Haux )
end

import Base: length
function length(e::ElecGradient)
    return length(e.psiks)
end

mutable struct ElecVars
    psiks::BlochWavefunc
    Hsub::Array{Matrix{ComplexF64},1}
    Haux_eigs::Array{Float64,2}
end

function ElecVars( Ham::Hamiltonian )
    return ElecVars( Ham, rand_BlochWavefunc(Ham) )
end

function ElecVars( Ham::Hamiltonian, psiks::BlochWavefunc )
    
    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin

    Hsub = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    Haux_eigs = zeros(Float64,Nstates,Nkspin) # the same as electrons.ebands
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        #
        Hsub[i] = zeros(ComplexF64,Nstates,Nstates)
        #
        Hsub[i][:] = psiks[i]' * op_H(Ham, psiks[i])
        #
        Haux_eigs[:,i] = eigvals(Hermitian(Hsub[i]))  # set Haux_eigs to eigenvalues of Hsub
    end

    return ElecVars(psiks, Hsub, Haux_eigs)
end

import Base: show
function show( io::IO, evars::ElecVars )
    Nkspin = length(evars.psiks)
    for i in 1:Nkspin
        println("Haux i = ", i)
        display(evars.Hsub[i]); println()
        println("Haux_eigs i = ", i)
        display( eigvals(Hermitian(evars.Hsub[i])) ); println()
    end
end
show( evars::ElecVars ) = show( stdout, evars )




