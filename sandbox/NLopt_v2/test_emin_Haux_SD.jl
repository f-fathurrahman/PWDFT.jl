using LinearAlgebra
using Printf
using PWDFT
using Random
import Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("calc_grad_Haux.jl")
include("grad_obj_function.jl")
include("create_Ham.jl")

function precond_grad!( Ham, g, Kg )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        Kg[ikspin] = Kprec( ik, Ham.pw, g[ikspin] )
    end # ikspin
    return
end


# steepest descent
function main_SD()

    Random.seed!(1234)

    #Ham = create_Ham_atom_Al_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()

    psiks = rand_BlochWavefunc( Ham )
    Rhoe = guess_rhoe_atomic( Ham )
    update!(Ham, Rhoe)
    _ = diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=30 )

    #psiks = Serialization.deserialize("psiks.data")

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

    Etot_old = obj_function!( Ham, psiks, Haux, skip_ortho=true )

    println("Etot_old = ", Etot_old)

    α_t = 1e-5

    for iter = 1:100
        
        grad_obj_function!( Ham, psiks, g, Haux, g_Haux )

        precond_grad!( Ham, g, Kg )
        Kg_Haux = g_Haux

        # update psiks
        psiks = psiks - α_t*g

        # update Haux
        for i = 1:Nkspin
            Haux[i] = Haux[i] - α_t*Kg_Haux[i]
            Haux[i] = 0.5*( Haux[i] + Haux[i]' ) # or use previous U_Haux
        end

        Etot = obj_function!( Ham, psiks, Haux )

        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        check_evals!( Ham, psiks, Haux, Hsub )

        Etot_old = Etot
    end

    Serialization.serialize( "psiks.data", psiks )

end

function check_evals!( Ham, psiks, Haux, Hsub )
    # Calculate Hsub for comparison with Haux
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    s = 0.0
    for ispin in 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin -1)*Nkpt
        Hsub[i] = psiks[i]' * ( Ham*psiks[i] )
        evalsHsub = eigvals(Hsub[i])
        evalsHaux = diag(Haux[i])
        s = s + abs(sum(evalsHaux) - sum(evalsHsub))
    end
    println("diff sum evals = ", s)
    return
end


main_SD()

