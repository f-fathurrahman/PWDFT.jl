using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")

function set_occupations!( Ham, kT )

    Ham.electrons.Focc, E_fermi = calc_Focc(
        Ham.electrons.Nelectrons,
        Ham.pw.gvecw.kpoints.wk,
        kT, Ham.electrons.ebands,
        Ham.electrons.Nspin )

    return E_fermi
end


function eval_L( Ham, psiks; kT=1e-3 )

    E_fermi = set_occupations!( Ham, kT )
    Entropy = calc_entropy(
        Ham.pw.gvecw.kpoints.wk,
        kT,
        Ham.electrons.ebands,
        E_fermi,
        Ham.electrons.Nspin )

    energies = calc_energies( Ham, psiks )
    energies.mTS = Entropy

    print_ebands( Ham.electrons, Ham.pw.gvecw.kpoints )

    return sum(energies)
end

# modify Ham.electrons.ebands
# rotate psiks
# return subspace Hamiltonian
# psiks should be orthonormalized first
function subspace_rotation!( Ham, psiks )

    Nkspin = length(psiks)
    Nstates = Ham.electrons.Nstates
    Hsub = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Hsub[i] = zeros(ComplexF64,Nstates,Nstates)
    end

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        Ham.electrons.ebands[:,i], evecs = eigen(Hr)
        psiks[i] = psiks[i]*evecs # also rotate
        Hsub[i] = psiks[i]' * ( op_H(Ham, psiks[i]) )
    end

    return Hsub
end

function test_main()
    Ham = create_Ham_atom_Pt_smearing()
    psiks = rand_BlochWavefunc(Ham)

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    Hsub = subspace_rotation!( Ham, psiks )
    
    Etot = eval_L(Ham, psiks)
    @printf("Etot = %18.10f\n", Etot)

    for i in 1:length(Hsub)
        println("Real part\n")
        display(real(Hsub[i]))
        println("")
        println("Imaginary part\n")
        display(imag(Hsub[i]))
        println("")
    end

    Etot = eval_L(Ham, psiks)
    @printf("Etot = %18.10f\n", Etot)

end

test_main()