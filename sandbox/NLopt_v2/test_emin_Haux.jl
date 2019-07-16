using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("calc_grad_Haux.jl")
include("grad_obj_function.jl")

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

function precond_grad!( Ham, g, Kg )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin-1)*Nkpt
            Kg[ikspin] = Kprec( ik, Ham.pw, g[ikspin] )
        end
    end
    return
end


# steepest descent
function main_SD()

    Random.seed!(1234)

    Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    
    if Ham.sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    end

    psiks = rand_BlochWavefunc( Ham )

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = length(psiks)

    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = rand( ComplexF64, Nstates, Nstates )
        Haux[i] = 0.5*( Haux[i] + Haux[i]' )
    end

    Hsub = copy(Haux)

    g = zeros_BlochWavefunc( Ham )
    Kg = copy(g)

    g_Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i = 1:Nkspin
        g_Haux[i] = zeros( ComplexF64, Nstates, Nstates )
    end
    Kg_Haux = copy(g_Haux)

    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true, rhoe_symm=rhoe_symmetrizer )
    #for i = 1:Nkspin
    #    println("\nikspin = ", i)
    #    println("\nreal part")
    #    print_matrix(real(Haux[i]))
    #    println("\nimag part")
    #    print_matrix(imag(Haux[i]))
    #end

    println("Etot_old = ", Etot_old)

    α_t = 1e-5

    for iter = 1:10
        
        grad_obj_function!( Ham, psiks, g, Haux, g_Haux, rhoe_symm=rhoe_symmetrizer )

        precond_grad!( Ham, g, Kg )
        Kg_Haux = 0.01*g_Haux

        psiks = psiks - α_t*g

        for i = 1:Nkspin
            #println("\nikspin = ", i)
            #println("\nreal part")
            #print_matrix(real(g_Haux[i]))
            #println("\nimag part")
            #print_matrix(imag(g_Haux[i]))
            Haux[i] = Haux[i] + α_t*Kg_Haux[i]
            Haux[i] = 0.5*( Haux[i] + Haux[i]' ) # or use previous U_Haux
        end

        Etot = obj_function!( Ham, psiks, Haux, rhoe_symm=rhoe_symmetrizer )

        #@printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        # Calculate Hsub (for comparison with Haux)
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            #
            Hsub[ikspin] = psiks[ikspin]' * ( Ham*psiks[ikspin] )

            evalsHsub = eigvals(Hsub[ikspin])
            evalsHaux = diag(Haux[ikspin])

            #println("diff Haux = ", sum(Hsub[ikspin] - Haux[ikspin]))
            println("diff sum evals = ", abs(sum(evalsHaux) - sum(evalsHsub)))

        end
        end


        Etot_old = Etot
    end

end
#main_SD()



function calc_beta_CG!( g, g_old, Kg, Kg_old, β )
    for i = 1:length(g)
        β[i] = real(sum(conj(g[i]-g_old[i]).*Kg[i]))/real(sum(conj(g_old[i]).*Kg_old[i]))
        if β[i] < 0.0
            β[i] = 0.0
        end
    end
    return
end


function calc_alpha_CG!( α_t, g, gt, d, α )
    for i = 1:length(g)
        denum = real(sum(conj(g[i]-gt[i]).*d[i]))
        if denum != 0.0
            α[i] = abs( α_t*real(sum(conj(g[i]).*d[i]))/denum )
        else
            α[i] = 0.0
        end
    end
    return
end


function main_CG()

    Random.seed!(1234)

    Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()

    if Ham.sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    end


    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = rand( ComplexF64, Nstates, Nstates )
        Haux[i] = 0.5*( Haux[i] + Haux[i]' )
    end

    Hsub = copy(Haux)




    psiks = rand_BlochWavefunc( Ham )

    # prepare guess wavefunc
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64,Npoints,Nspin)
    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )
    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
    end
    #
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )



    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )


    g  = zeros_BlochWavefunc( Ham )
    gt = copy(g) 
    Kg = copy(g)
    g_old = copy(g)
    Kg_old = copy(g)
    d = copy(g)
    d_old = copy(g)
    psic = copy(g)

    g_Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i = 1:Nkspin
        g_Haux[i] = zeros( ComplexF64, Nstates, Nstates )
    end
    Kg_Haux     = copy(g_Haux)
    gt_Haux     = copy(g_Haux) 
    g_Haux_old  = copy(g_Haux)
    Kg_Haux_old = copy(g_Haux)
    d_Haux      = copy(g_Haux)
    d_Haux_old  = copy(g_Haux)
    Hauxc       = copy(g_Haux)

    β = zeros(length(g))
    α = zeros(length(g))

    β_Haux = zeros(length(g_Haux))
    α_Haux = zeros(length(g_Haux))



    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    α_t = 3e-5
    etot_conv_thr = 1e-6

    Nconverges = 0

    for iter = 1:50
        
        grad_obj_function!( Ham, psiks, g, Haux, g_Haux, rhoe_symm=rhoe_symmetrizer )
        precond_grad!( Ham, g, Kg )

        Kg_Haux = 0.05*g_Haux  # scalar preconditioner

        if iter > 1
            calc_beta_CG!( g, g_old, Kg, Kg_old, β )
            calc_beta_CG!( g_Haux, g_Haux_old, Kg_Haux, Kg_Haux_old, β_Haux )
        end

        d = -Kg + β .* d_old

        d_Haux = Kg_Haux + β_Haux .* d_Haux_old

        psic = psiks + α_t*d  # trial wavefunc

        Hauxc = Haux + α_t*d_Haux
        for i in 1:Nkspin
            Hauxc[i] = 0.5*( Hauxc[i] + Hauxc[i]' )
        end

        # line minimization
        grad_obj_function!( Ham, psic, gt, Hauxc, gt_Haux, rhoe_symm=rhoe_symmetrizer )
        
        calc_alpha_CG!( α_t, g, gt, d, α )
        psiks = psiks + α .* d

        calc_alpha_CG!( α_t, g_Haux, gt_Haux, d_Haux, α_Haux )
        for i in 1:Nkspin
            Haux[i] = Haux[i] + α_Haux[i] .* d_Haux[i]
            Haux[i] = 0.5( Haux[i] + Haux[i]' )
        end

        Etot = obj_function!( Ham, psiks, Haux, rhoe_symm=rhoe_symmetrizer )

        diffE = abs(Etot_old - Etot)
        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        # Calculate Hsub (for comparison with Haux)
        #=for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            #
            Hsub[ikspin] = psiks[ikspin]' * ( Ham*psiks[ikspin] )

            evalsHsub = eigvals(Hsub[ikspin])
            evalsHaux = diag(Haux[ikspin])

            #println("diff Haux = ", sum(Hsub[ikspin] - Haux[ikspin]))
            println("diff sum evals = ", abs(sum(evalsHaux) - sum(evalsHsub)))

        end
        end=#

        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nEmin_PCG is converged in iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        Etot_old = Etot
        g_old = copy(g)
        Kg_old = copy(Kg)
        d_old = copy(d)

        g_Haux_old = copy(g_Haux)
        Kg_Haux_old = copy(Kg_Haux)
        d_Haux_old = copy(d_Haux)


    end

end

@time main_CG()

