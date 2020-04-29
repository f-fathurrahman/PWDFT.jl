using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("MinimizeParams.jl")
include("test_ElecVars.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("linmin_quad.jl")
include("setup_guess_wavefunc.jl")

function main()

    Random.seed!(1234)

    kT = 0.01
    #Ham = create_Ham_atom_Si_smearing()
    Ham = create_Ham_atom_Al_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()

    #psiks = rand_BlochWavefunc(Ham)
    #setup_guess_wavefunc!( Ham, psiks, startingrhoe=:gaussian, skip_initial_diag=false )
    #evars = ElecVars(Ham, psiks)
    
    evars = ElecVars(Ham)
    
    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    Nkspin = Nkpt*Nspin

    subrot = SubspaceRotations(Nkspin, Nstates)

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    Etot_old = compute!( Ham, evars, g, Kg, kT, subrot )
    println("Etot_old = ", Etot_old)
    psiks_orig = deepcopy(evars.psiks)
    Hsub_orig = deepcopy(evars.Hsub)
    println("Hsub orig")
    print_vec_mat(Hsub_orig)
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)

    d = deepcopy(g)
    for i in 1:Nkspin
        d.psiks[i] = -Kg.psiks[i]
        d.Haux[i]  = -Kg.Haux[i] 
    end
    constrain_search_dir!( d, evars )

    #α, α_Haux = linmin_grad_v2!( Ham, evars, g, d, kT, subrot )
    #@printf("α = %f, α_Haux = %f\n", α, α_Haux)
    α = 0.0
    α_Haux = 1.0
    do_step!( α, α_Haux, evars, d, subrot )

    psiks_now = deepcopy(evars.psiks)
    Hsub_now = deepcopy(evars.Hsub)
    println("Hsub now")
    print_vec_mat(Hsub_now)
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)

    Etot = compute!( Ham, evars, g, Kg, kT, subrot )
    println("Etot = ", Etot)
    diffE = Etot - Etot_old
    if diffE > 0
        println("Energy is not reducing")
    else
        println()
    end
    print_vec_mat(evars.Hsub)
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)

end

main()