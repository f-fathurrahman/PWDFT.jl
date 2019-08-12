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
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()

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
    @printf("Etot_old_v1 = %18.10f\n", Etot_old_v1)

    Etot_old_v1 = obj_function_v1!( Ham, psiks, skip_ortho=true )
    @printf("Etot_old_v1 = %18.10f\n", Etot_old_v1)

    Ham.electrons.ebands = rotate_subspace( Ham, psiks, skip_ortho=true )
    for i in 1:Nkspin
        Haux[i] = diagm(0 => Ham.electrons.ebands[:,i])
    end
    
    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )
    @printf("Etot_old    = %18.10f\n", Etot_old)

    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )
    @printf("Etot_old    = %18.10f\n", Etot_old)

    g = copy(psiks)
    g_Haux = copy(Haux)

    for itry = 1:5

        println("\nitry = ", itry)

        grad_obj_function!( Ham, psiks, g, Haux, g_Haux )
        s1 = sum(g)
        ss1 = sum(g_Haux)
        @printf("sum g      = [%18.10f,%18.10f]\n", real(s1), imag(s1))
        @printf("sum g_Haux = [%18.10f,%18.10f]\n", real(ss1), imag(ss1))

        grad_obj_function!( Ham, psiks, g, Haux, g_Haux )
        s2 = sum(g)
        ss2 = sum(g_Haux)
        @printf("sum g      = [%18.10f,%18.10f]\n", real(s2), imag(s2))
        @printf("sum g_Haux = [%18.10f,%18.10f]\n", real(ss2), imag(ss2))

        @printf("diff g      = [%18.10f,%18.10f]\n", real(s1 - s2), imag(s1 - s2))
        @printf("diff g_Haux = [%18.10f,%18.10f]\n", real(ss1 - ss2), imag(ss1 - ss2))
    end
end


import Base: sum
function sum( f::Array{Array{ComplexF64,2},1} )
    s = 0.0 + im*0.0
    for i in 1:length(f)
        s = s + sum(f[i])
    end
    return s
end

function rotate_subspace( Ham, psiks_; skip_ortho=false )
    psiks = copy(psiks_)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    ebands = zeros(Float64, Nstates, Nkspin)

    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        ebands[:,i], evecs = eigen(Hr)
    end

    return ebands
end

test_main()