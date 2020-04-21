mutable struct ElecVars
    psiks::BlochWavefunc
    Hsub::Array{Matrix{ComplexF64},1}
    Haux_eigs::Array{Float64,2}
end

mutable struct ElecGradient
    psiks::BlochWavefunc
    Haux::Array{Matrix{ComplexF64},1}
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

function ElecVars(Ham::Hamiltonian)

    psiks = rand_BlochWavefunc(Ham)
    
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





function update_occ!( Ham, evars::ElecVars, kT )
    
    Nspin = Ham.electrons.Nspin
    Nelectrons = Ham.electrons.Nelectrons
    wk = Ham.pw.gvecw.kpoints.wk

    Ham.electrons.ebands = copy(evars.Haux_eigs)
    evals = evars.Haux_eigs

    Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
    mTS = calc_entropy( wk, kT, evals, E_fermi, Nspin )
    
    Ham.electrons.Focc = copy(Focc)

    return E_fermi, mTS

end

"""
Return total energies
"""
function calc_energies_grad!( Ham::Hamiltonian, evars::ElecVars, egrad::ElecGradient, kT::Float64 )

    psiks = evars.psiks
    
    update_occ!( Ham, evars::ElecVars, kT )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    
    Ham.energies = calc_energies( Ham, psiks )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        calc_grad!( Ham, psiks[i], g[i] )
        Kprec!( ik, Ham.pw, g[i], Kg[i] )
    end

    return sum( Ham.energies )
end