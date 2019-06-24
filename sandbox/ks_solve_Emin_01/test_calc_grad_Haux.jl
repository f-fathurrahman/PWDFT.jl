using Random
using LinearAlgebra
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include("calc_grad_Haux.jl")

function create_Ham_Pt_fcc_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 30.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[3,3,3], extra_states=4 )
end

function create_Ham_atom_Pt_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 30.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_Al_fcc_smearing()
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(7.6525970200) )
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[3,3,3], extra_states=4 )

end


function main()
    Random.seed!(1234)
    
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    Ham = create_Ham_atom_Pt_smearing()

    psiks = rand_BlochWavefunc(Ham)

    # random Haux
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nelectrons = Ham.electrons.Nelectrons
    wk = Ham.pw.gvecw.kpoints.wk

    Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    eta = zeros(Float64,Nstates,Nkspin)
    for ikspin = 1:Nkspin
        Haux[ikspin] = rand(ComplexF64, Nstates, Nstates)
        Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]')
        eta[:,ikspin] = eigvals(Haux[ikspin])
    end

    g_psi = zeros_BlochWavefunc( Ham )

    g_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Kg_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Haux_c = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ikspin = 1:Nkspin
        g_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        Haux_c[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux_old[ikspin] = zeros(ComplexF64,Nstates,Nstates)        
    end

    β_Haux = zeros(Float64,Nkspin)

    kT = 0.001

    Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
    Ham.electrons.Focc[:,:] = Focc

    α_t = 3e-5

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt

        g_psi[ikspin], g_Haux[ikspin] = calc_grad_Haux( Ham, psiks[ikspin], eta[:,ikspin], kT )

        Kg_Haux[ikspin] = 0.1*g_Haux[ikspin]

        #println("Real Matrix:")
        #print_matrix(real(g_Haux[ikspin]))
        #println("Imag Matrix:")
        #print_matrix(imag(g_Haux[ikspin]))

        #β[ikspin] =
        #real(sum(conj(g_Haux[ikspin]-g_Haux_old[ikspin]).*Kg_Haux[ikspin]))/
        #real(sum(conj(g_Haux_old[ikspin]).*Kg_Haux_old[ikspin]))

        d_Haux[ikspin] = -Kg_Haux[ikspin] + β_Haux[ikspin] * d_Haux_old[ikspin]
        
        Haux_c[ikspin] = Haux[ikspin] + α_t*d_Haux[ikspin]

        println("Real Haux")
        print_matrix(real(Haux_c[ikspin]))
        println("Imag Haux")
        print_matrix(imag(Haux_c[ikspin]))

        Haux_c[ikspin] = 0.5*(Haux_c[ikspin] + Haux_c[ikspin]') # make the eigenvalues real
        eta[:,ikspin] = eigvals(Haux_c[ikspin])
        println("Eigenvalues")
        for i = 1:Nstates
            println( eta[i,ikspin] )
        end

    end
    end

end

main()