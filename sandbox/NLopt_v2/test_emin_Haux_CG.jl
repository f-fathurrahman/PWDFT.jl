using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("calc_grad_Haux.jl")
include("grad_obj_function.jl")
include("create_Ham.jl")

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
        #if abs(denum) <= 1e-6
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

    #Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_Pt_fcc_smearing()

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
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )

    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )
    println("Etot_old = ", Etot_old)

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
        
        grad_obj_function!( Ham, psiks, g, Haux, g_Haux )
        precond_grad!( Ham, g, Kg )

        Kg_Haux = 0.01*g_Haux  # scalar preconditioner

        if iter > 1
            calc_beta_CG!( g, g_old, Kg, Kg_old, β )
            calc_beta_CG!( g_Haux, g_Haux_old, Kg_Haux, Kg_Haux_old, β_Haux )
        end

        d = -Kg + β .* d_old
        # line minimization
        psic = psiks + α_t*d  # trial wavefunc

        #d_Haux = Kg_Haux + β_Haux .* d_Haux_old
        #Hauxc = Haux + α_t*d_Haux
        #for i in 1:Nkspin
        #    Hauxc[i] = 0.5*( Hauxc[i] + Hauxc[i]' )
        #end

        grad_obj_function!( Ham, psic, gt, Hauxc, gt_Haux )       
        
        calc_alpha_CG!( α_t, g, gt, d, α )
        psiks = psiks + α .* d

        #calc_alpha_CG!( α_t, g_Haux, gt_Haux, d_Haux, α_Haux )
        #for i in 1:Nkspin
        #    println("α_Haux = ", α_Haux[i])
        #    Haux[i] = Haux[i] + α_Haux[i] .* d_Haux[i]
        #    Haux[i] = 0.5( Haux[i] + Haux[i]' )
        #end

        Etot = obj_function!( Ham, psiks, Haux )

        diffE = abs(Etot_old - Etot)
        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

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

