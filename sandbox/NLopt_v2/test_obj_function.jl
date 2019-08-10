using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("obj_function_v1.jl")
include("calc_grad_Haux.jl")
include("grad_obj_function.jl")
include("create_Ham.jl")

function test_main()

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

    psiks_old = copy(psiks)

    Etot_old_v1 = obj_function_v1!( Ham, psiks, skip_ortho=true )
    
    #update!( Ham, Rhoe_old )
    #Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )
    
    #println("Etot_old    = ", Etot_old)
    println("Etot_old_v1 = ", Etot_old_v1)

    #update!( Ham, Rhoe_old )
    Etot_old_v1 = obj_function_v1!( Ham, psiks, skip_ortho=true )

    #update!( Ham, Rhoe_old )
    #Etot_old = obj_function!( Ham, psiks_old, Haux, skip_ortho=true )
    
    #println("Etot_old    = ", Etot_old)
    println("Etot_old_v1 = ", Etot_old_v1)

    #for i in 1:Nkspin
    #    Haux[i] = diagm(0 => Ham.electrons.ebands[:,i])
    #end
    #Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )
    #println("Etot_old    = ", Etot_old)

end

test_main()