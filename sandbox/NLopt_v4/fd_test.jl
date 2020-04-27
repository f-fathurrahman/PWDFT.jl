using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("test_ElecVars.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("setup_guess_wavefunc.jl")
include("KS_solve_Emin_SD_Haux.jl")
include("KS_solve_Emin_PCG_Haux_v1.jl")
include("KS_solve_Emin_PCG_Haux_v2.jl")

function main()
    Random.seed!(1234)

    kT = 0.01
    #Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    #println(Ham)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates
    Nkspin = Nkpt*Nspin

    rotPrev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        rotPrev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end

    evars = ElecVars( Ham )
    g = ElecGradient( Ham )
    Kg = ElecGradient( Ham )

    Ham.energies.NN = calc_E_NN( Ham.atoms )

    E0 = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )

    println("E0 = ", E0)

    #d = -Kg
    d = deepcopy(g)

    constrain_search_dir!( d, evars )

    println("dot(evars.psiks,g.psiks) = ", dot(evars.psiks, g.psiks))
    Haux = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    g_Haux = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = diagm( 0 => evars.Haux_eigs[:,i] )
        g_Haux[i]= diagm( 0 => diag(g.Haux[i]) )
    end
    println("dot(evars.Haux,g.Haux) = ", dot(Haux, g.Haux))
    println("dot(evars.Haux,g.Haux) = ", dot(Haux, g_Haux))

    exit()
    #println("dot(psiks,psiks) = ", dot(evars.psiks, g.psiks)) # dot for Haux?
    #println("dot(d,d) = ", dot(d,d))

    scaled_d = deepcopy(d)

    psic = zeros_BlochWavefunc( Ham )
    for α in 10.0 .^ range(1,stop=-10,step=-1)
        
        scaled_d = ElecGradient(α*d.psiks, α*d.Haux)
        dE = dot_ElecGradient(g, scaled_d)
        
        @printf("α = %e, dE = %18.10e\n", α, dE)
        #
        evarsc = deepcopy(evars)
        rotPrevc = deepcopy(rotPrev)
        rotPrevCc = deepcopy(rotPrevC)
        rotPrevCinvc = deepcopy(rotPrevCinv)
        gt = deepcopy(g)
        Kgt = deepcopy(Kg)

        do_step!( α, evarsc, d, rotPrevc, rotPrevCc, rotPrevCinvc )
        
        Etot = compute!( Ham, evarsc, gt, Kgt, kT, rotPrevCinvc, rotPrevc ) # 
        
        ratio = (Etot - E0)/dE
        @printf("α = %e, ratio = %18.10e, diff = %e\n", α, ratio, abs(1-ratio))
    end

    println("Pass here")

end

main()