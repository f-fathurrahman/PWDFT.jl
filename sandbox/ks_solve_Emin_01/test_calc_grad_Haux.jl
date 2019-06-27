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
    ecutwfc = 15.0
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

function create_Ham_atom_Al_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end


function main()
    Random.seed!(1234)
    
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_atom_Al_smearing()    

    pw = Ham.pw

    # random Haux
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nelectrons = Ham.electrons.Nelectrons
    wk = Ham.pw.gvecw.kpoints.wk
    Npoints = prod(Ham.pw.Ns)

    Rhoe = zeros(Float64,Npoints,Nspin)

    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )

    update!(Ham, Rhoe)

    psiks = rand_BlochWavefunc( Ham )

    Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    eta = zeros(Float64,Nstates,Nkspin)
    for ikspin = 1:Nkspin
        Haux[ikspin] = rand(ComplexF64, Nstates, Nstates)
        Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]')
        eta[:,ikspin] = eigvals(Haux[ikspin])
    end

    g_psi = zeros_BlochWavefunc( Ham )
    g_psi_old = zeros_BlochWavefunc( Ham )
    Δ_psi = zeros_BlochWavefunc( Ham )
    Δ_psi_old = zeros_BlochWavefunc( Ham )

    Hsub = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    Δ_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    g_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Kg_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Haux_c = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ikspin = 1:Nkspin

        Hsub[ikspin] = zeros(ComplexF64,Nstates,Nstates)

        Δ_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        g_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        Haux_c[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux_old[ikspin] = zeros(ComplexF64,Nstates,Nstates)        
    end

    β_Haux = zeros(Float64,Nkspin)

    λ_psi = zeros(Float64,Nkspin)
    λ_Haux = zeros(Float64,Nkspin)

    kT = 0.001

    Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
    Ham.electrons.Focc[:,:] = Focc

    α_t = 3e-4
    Etot_old = 1.0

    for iter = 1:50

        for ispin = 1:Nspin
        for ik = 1:Nkpt
            
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
    
            g_psi[ikspin], g_Haux[ikspin], Δ_psi[ikspin], Δ_Haux[ikspin] =
            calc_grad_Haux( Ham, psiks[ikspin], eta[:,ikspin], kT )
    
            Δ_psi[ikspin] = -Kprec( Ham.ik, pw, Δ_psi[ikspin] )
            Δ_Haux[ikspin] = 0.1*Δ_Haux[ikspin]

        #if iter > 1
        #    λ_psi[ikspin] =
        #    real(sum(conj(g_psi[ikspin]).*Δ_psi[ikspin]))/
        #    real(sum(conj(g_psi_old[ikspin])).*Δ_psi_old[ikspin] )
        #    #
        #    λ_Haux[ikspin] =
        #    real(sum(conj(g_Haux[ikspin]).*Δ_Haux[ikspin]))/
        #    real(sum(conj(g_Haux_old[ikspin])).*Δ_Haux_old[ikspin] )            
        #end
#
#        #d_psi[ikspin] = Δ_psi[ikspin] + λ_psi[ikspin] * d_psi_old[ikspin]
        #d_Haux[ikspin] = Δ_Haux[ikspin] + λ_Haux[ikspin] * d_Haux_old[ikspin]

        #Haux_c[ikspin] = 0.5*(Haux_c[ikspin] + Haux_c[ikspin]') # make the eigenvalues real
        #eta[:,ikspin] = eigvals(Haux_c[ikspin])

            psiks[ikspin] = psiks[ikspin] + 1e-5*Δ_psi[ikspin]
            ortho_sqrt(psiks[ikspin])

            Haux[ikspin] = Haux[ikspin] + 1e-3*Δ_Haux[ikspin]
            Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]') # make the eigenvalues real
            eta[:,ikspin] = eigvals(Haux[ikspin])
        
            #println("Eigenvalues")
            #for i = 1:Nstates
            #    println( eta[i,ikspin] )
            #end

            #Δ_psi_old = copy(Δ_psi)
            #g_psi_old = copy(g_psi)

        end
        end
    
        Ham.electrons.Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
        Entropy = calc_entropy( wk, kT, eta, E_fermi, Nspin )

        calc_rhoe!( Ham, psiks, Rhoe )
        update!(Ham, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Ham.energies.mTS = Entropy
        Etot = sum(Ham.energies)

        # Calculate Hsub (for comparison with Haux)
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            #
            Hsub[ikspin] = psiks[ikspin]' * ( Ham*psiks[ikspin] )

            println("diff Haux = ", real(sum(Hsub[ikspin] - Haux[ikspin])))
        end
        end



        @printf("%18.10f %18.10e\n", Etot, Etot_old-Etot)

        Etot_old = Etot

    end

end

main()