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

include("unfold_BlochWavefuncGamma.jl")

function test_01()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )

    psis = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis)

    Nstates = size(psis.data[1],2)
    ispin = 1
    for ist in 1:Nstates
        println()
        for igw in 1:3
            c = psis.data[ispin][igw,ist]
            @printf("%3d %3d [%18.10f,%18.10f]\n", igw, ist, c.re, c.im)
        end
    end

    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    Rhoe = calc_rhoe(Ham, psis)

    Rhoe_ = calc_rhoe(Ham_, psiks)

    update!(Ham, Rhoe)
    update!(Ham_, Rhoe_)

    println("V Ps loc comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Ps_loc[ip], Ham_.potentials.Ps_loc[ip])
    end

    println("V Hartree comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Hartree[ip], Ham_.potentials.Hartree[ip])
    end

    println("V XC comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.XC[ip,1], Ham_.potentials.XC[ip,1])
    end

    println("V Total comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Total[ip,1], Ham_.potentials.Total[ip,1])
    end

    println()
    Kpsis = op_K(Ham, psis)
    s1 = sum(Kpsis.data[1])
    println("sum Kpsis  = ", s1 + conj(s1))

    Kpsis_ = op_K(Ham_, psiks)
    s2 = sum(Kpsis_[1])
    println("sum Kpsis_ = ", s2)

    println()

    V_loc_psis = op_V_loc(Ham, psis)
    s1 = sum(V_loc_psis.data[1])
    s1 = s1 + conj(s1)
    println("sum V_loc_psis  = ", s1)

    V_loc_psis_ = op_V_loc(Ham_, psiks)
    s2 = sum(V_loc_psis_[1])
    println("sum V_loc_psis_ = ", s2)

    println("s1 - s2    = ", s1 - s2)
    ss = 0.0 + 0.0*im
    for ist in 1:size(psis.data[1],2)
        ss = ss + V_loc_psis.data[1][1,ist]
    end
    println("V_loc_psis = ", ss)

    ist = 1
    for igw in 1:5
        c1 = V_loc_psis.data[1][igw,ist]
        c2 = V_loc_psis_[1][igw,ist]
        @printf("%3d [%18.10f,%18.10f] [%18.10f,%18.10f]\n", igw, c1.re, c1.im, c2.re, c2.im)
    end

    println(dot(psis,psis))
    println(dot(psiks,psiks))

    println("dot V_loc_psis")
    println(dot(V_loc_psis,V_loc_psis))
    println(dot(V_loc_psis_,V_loc_psis_))
    println(dot_BlochWavefuncGamma(V_loc_psis,V_loc_psis))

    V_Ps_nloc_psis = op_V_Ps_nloc(Ham, psis)
    V_Ps_nloc_psiks = op_V_Ps_nloc(Ham_, psiks)

    println("dot V_Ps_nloc_psis")
    println( dot(V_Ps_nloc_psis,V_Ps_nloc_psis) )
    println( dot(V_Ps_nloc_psiks,V_Ps_nloc_psiks) )
    println( dot_BlochWavefuncGamma(V_Ps_nloc_psis,V_Ps_nloc_psis) )

    ist = 1
    for igw in 1:5
        c1 = V_Ps_nloc_psis.data[1][igw,ist]
        c2 = V_Ps_nloc_psiks[1][igw,ist]
        @printf("%3d [%18.10f,%18.10f] [%18.10f,%18.10f]\n", igw, c1.re, c1.im, c2.re, c2.im)
    end

    Hpsis = op_H(Ham, psis)
    Hpsiks = op_H(Ham_, psiks)

    println("dot Hpsis")
    println( "dot      = ", dot(Hpsis, Hpsis) )
    println( "dot_orig = ", dot_orig(Hpsis, Hpsis) )
    println( dot(Hpsiks,Hpsiks) )
    println( dot_BlochWavefuncGamma(Hpsis, Hpsis) )

    Hsub = psis.data[1]' * Hpsis.data[1]
    Hsub = Hsub + conj(Hsub)

    Hsub_ = psiks[1]' * Hpsiks[1]

    println()
    println("Hsub real part  = "); display(real(Hsub)); println()
    println("Hsub_ real part = "); display(real(Hsub_)); println()
    
    println("Hsub imag part  = "); display(imag(Hsub)); println()
    println("Hsub_ imag part = "); display(imag(Hsub_)); println()



end

test_01()