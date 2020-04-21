using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("test_ElecVars.jl")
include("calc_energies_grad.jl")
include("emin_smearing.jl")

function main()

    Random.seed!(1234)

    kT = 0.01
    #Ham = create_Ham_atom_Al_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    println(Ham)

    #test_ElecVars(Ham)

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    evars = ElecVars(Ham)
    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)

    rotPrev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        rotPrev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end


    Ham.energies.NN = calc_E_NN(Ham.atoms)
    
    #Etot = calc_energies_grad!( Ham, evars, g, Kg, kT )
    Etot = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )

    println(Ham.energies)
    println("Etot = ", Etot)

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, evars )


    println("Pass here")
end

main()