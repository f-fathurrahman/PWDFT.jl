using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
include("BlochWavefuncGamma.jl")
include("calc_rhoe_gamma.jl")
include("Poisson_solve_gamma.jl")
include("op_K_gamma.jl")
include("op_V_loc_gamma.jl")
include("op_V_Ps_nloc_gamma.jl")
include("op_H_gamma.jl")
include("calc_energies_gamma.jl")
include("calc_grad_gamma.jl")

include("setup_guess_wavefunc.jl")

include("KS_solve_Emin_PCG_dot.jl")
include("calc_energies_grad.jl")
include("linmin_grad.jl")

include("KS_solve_Emin_PCG_dot_gamma.jl")
include("calc_energies_grad_gamma.jl")
include("linmin_grad_gamma.jl")

include("unfold_BlochWavefuncGamma.jl")

function test_01()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    #
    # Initialize Hamiltonian
    #
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    #
    # Wavefunctions
    #
    psis = randn_BlochWavefuncGamma(Ham)
    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    # Should be the same for non-gamma and gamma
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_old = zeros(Float64,Npoints,Nspin)

    Rhoe_ = zeros(Float64,Npoints,Nspin)


    g     = zeros_BlochWavefuncGamma(Ham)
    Kg    = zeros_BlochWavefuncGamma(Ham)
    gPrev = zeros_BlochWavefuncGamma(Ham)
    Hsub = Vector{Matrix{ComplexF64}}(undef,Nspin)
    for ispin in 1:Nspin
        Hsub[ispin] = zeros(ComplexF64,Nstates,Nstates)
    end

    g_     = zeros_BlochWavefunc(Ham_)
    Kg_    = zeros_BlochWavefunc(Ham_)
    gPrev_ = zeros_BlochWavefunc(Ham_)


    startingrhoe = :random
    skip_initial_diag = true

    setup_guess_wavefunc!( Ham, psis, startingrhoe, skip_initial_diag )
    setup_guess_wavefunc!( Ham_, psiks, startingrhoe, skip_initial_diag )

    # calculate E_NN
    Ham.energies.NN  = calc_E_NN( Ham.atoms )
    Ham_.energies.NN = calc_E_NN( Ham_.atoms )

    println()
    println("Initial dot(psis,psis)   = ", 2*real(dot_BlochWavefuncGamma(psis,psis)))
    println("Initial dot(psiks,psiks) = ", dot_BlochWavefunc(psiks,psiks))

    # No need to orthonormalize
    Etot = calc_energies_grad!( Ham, psis, g, Kg, Hsub )
    Etot_ = calc_energies_grad!( Ham_, psiks, g_, Kg_ )

    println()
    println("Initial dot(g,g)   = ", 2*real(dot_BlochWavefuncGamma(g,g)))
    println("Initial dot(g_,g_) = ", dot_BlochWavefunc(g_,g_))

    println("Etot  = ", Etot)
    println("Etot_ = ", Etot_)

    println("\nSome g")
    println(g.data[1][1,1])
    println(g_[1][1,1])

    #println("Initial dot(psis,psis) = ", 2*real(dot_BlochWavefuncGamma(psis,psis)))    
    #println("Initial dot(g,g) = ", 2*real(dot_BlochWavefuncGamma(g,g)))
    #println("Initial dot(Kg,Kg) = ", 2*real(dot_BlochWavefuncGamma(Kg,Kg)))

    d = deepcopy(Kg)
    d_ = deepcopy(Kg_)
    
    #println("Before constrain_search_dir dot(d,d)   = ", 2*real(dot_BlochWavefuncGamma(d,d)))
    #println("Before constrain_search_dir dot(d_,d_) = ", dot_BlochWavefunc(d_,d_))
    #println("sum d = ", sum(d.data[1]))

    # Constrain
    constrain_search_dir!( d, psis )
    constrain_search_dir!( d_, psiks )

    #println("After constrain_search_dir dot(d,d)   = ", 2*real(dot_BlochWavefuncGamma(d,d)))
    #println("After constrain_search_dir dot(d_,d_) = ", dot_BlochWavefunc(d_,d_))
    #println("sum d = ", sum(d.data[1]))

    α = 0.0
    α_ = 0.0

    β  = 0.0
    β_ = 0.0

    psic = zeros_BlochWavefuncGamma(Ham)
    gt   = zeros_BlochWavefuncGamma(Ham)

    psic_ = zeros_BlochWavefunc(Ham_)
    gt_   = zeros_BlochWavefunc(Ham_)

    αt = 1e-5

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)
    dVol2 = Ham_.pw.CellVolume/prod(Ham_.pw.Ns)

    Etot_old = Etot
    Etot_old_ = Etot_

    for iter in 1:2

        println()
        @printf("================\n")
        @printf("Begin iter = %3d\n", iter)
        @printf("================\n")

        # Update search direction
        for i in 1:Nspin
            d.data[i] = -Kg.data[i] + β*d.data[i]
            d_[i] = -Kg_[i] + β_*d_[i]
        end

        println("Before constrain_search_dir dot(d,d)   = ", 2*real(dot_BlochWavefuncGamma(d,d)))
        println("Before constrain_search_dir dot(d_,d_) = ", dot_BlochWavefunc(d_,d_))
        println("sum d = ", sum(d.data[1]))

        # Constrain
        constrain_search_dir!( d, psis )
        constrain_search_dir!( d_, psiks )

        println("After constrain_search_dir dot(d,d)   = ", 2*real(dot_BlochWavefuncGamma(d,d)))
        println("After constrain_search_dir dot(d_,d_) = ", dot_BlochWavefunc(d_,d_))
        println("sum d = ", sum(d.data[1]))

        for i in 1:Nspin
            #
            psic.data[i][:] = psis.data[i] + αt*d.data[i]
            ortho_GS_gamma!(psic.data[i])
            #
            psic_[i] = psiks[i] + αt*d_[i]
            ortho_gram_schmidt!(psic_[i])
        end

        println("dot(psic)  = ", 2*dot_BlochWavefuncGamma(psic,psic))
        println("dot(psic_) = ", dot_BlochWavefunc(psic_,psic_))

        #
        # Calculate Rhoe
        #
        calc_rhoe!( Ham, psic, Rhoe )
        update!( Ham, Rhoe )

        calc_rhoe!( Ham_, psic_, Rhoe_ )
        update!( Ham_, Rhoe_ )
    
        println("integ Rhoe  = ", sum(Rhoe)*dVol)
        println("integ Rhoe_ = ", sum(Rhoe_)*dVol)

        for i in 1:Nspin
            #
            Ham.ispin = i
            calc_grad!( Ham, psic.data[i], gt.data[i] )
            #
            Ham_.ispin = i
            calc_grad!( Ham_, psic_[i], gt_[i] )
        end

        println("dot(gt)  = ", 2*dot_BlochWavefuncGamma(gt,gt))
        println("dot(gt_) = ", dot_BlochWavefunc(gt_,gt_))

        println("\nSome psic")
        println(psic.data[1][1,1])
        println(psic_[1][1,1])

        println("\nSome gt")
        println(gt.data[1][1,1])
        println(gt_[1][1,1])

        denum  = 2.0*real( dot_BlochWavefuncGamma(g - gt, d) )
        #c1 = dot_BlochWavefuncGamma(g, d)
        #c2 = dot_BlochWavefuncGamma(gt, d)
        #denum = 2*real(c1 - c2)
        
        denum_ = 2.0*real( dot(g_ .- gt_, d_) )

        println()
        println("denum  = ", denum)
        println("denum_ = ", denum_)
        println("diff = ", denum - denum_)

        if denum != 0.0
            α = abs( αt * 2.0*real( dot_BlochWavefuncGamma(g, d) )/denum )
        else
            α = 0.0
        end

        if denum_ != 0.0
            α_ = abs( αt * 2.0*real( dot(g_, d_) )/denum_ )
        else
            α_ = 0.0
        end

        println()
        println("α  = ", α)
        println("α_ = ", α_)
        println("diff = ", α - α_)

        # Update wavefunction
        for i in 1:Nspin
            #
            psis.data[i][:] = psis.data[i] + α*d.data[i]
            ortho_GS_gamma!( psis.data[i] )
            #
            psiks[i] = psiks[i] + α_*d_[i]
            ortho_gram_schmidt!( psiks[i] )
        end

        println("psis ortho check")

        ortho_check(psis)
        #
        println("dot(psis)  = ", 2*dot_BlochWavefuncGamma(psis,psis))
        println("dot(psiks) = ", dot_BlochWavefunc(psiks,psiks))

        println()
        println("Some psis and psiks")
        println(psis.data[1][1,1])
        println(psiks[1][1,1])

        Etot  = calc_energies_grad!( Ham, psis, g, Kg, Hsub )
        Etot_ = calc_energies_grad!( Ham_, psiks, g_, Kg_, )

        println()
        println("dot(g,g)   = ", 2*real(dot_BlochWavefuncGamma(g,g)))
        println("dot(g_,g_) = ", dot_BlochWavefunc(g_,g_))

        println("\nSome g")
        println(g.data[1][1,1])
        println(g_[1][1,1])

        dEtot  = Etot_old  - Etot
        dEtot_ = Etot_old_ - Etot_

        @printf("Iter = %3d, Etot  = %18.10f dEtot  = %18.10e\n", iter, Etot, dEtot)
        @printf("Iter = %3d, Etot_ = %18.10f dEtot_ = %18.10e\n", iter, Etot_, dEtot_)

        if dEtot < 0.0
            println("***WARNING dEtot is not decreasing")
        end

        if dEtot_ < 0.0
            println("***WARNING dEtot_ is not decreasing")
        end

        Etot_old  = Etot
        Etot_old_ = Etot_

    end

end

test_01()
