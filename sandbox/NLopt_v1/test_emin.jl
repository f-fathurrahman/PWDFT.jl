using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("grad_obj_function.jl")

function create_Ham_H2()
    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function create_Ham_H_atom()
    atoms = Atoms(xyz_string=
        """
        1

        H      0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end


# steepest descent
function main_SD()

    Random.seed!(1234)

    #Ham = create_Ham_H2()
    Ham = create_Ham_H_atom()

    psiks = rand_BlochWavefunc( Ham )
    Etot_old = obj_function!( Ham, psiks, skip_ortho=true )

    g = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )

    α_t = 1e-5

    for iter = 1:50
        grad_obj_function!( Ham, psiks, g, Kg )
        psiks = psiks - α_t*Kg
        Etot = obj_function!( Ham, psiks )
        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        Etot_old = Etot
    end

end
#main_SD()

function calc_beta_CG!( g, g_old, Kg, Kg_old, β )
    for i = 1:length(g)
        β[i] = real(sum(conj(g[i]-g_old[i]).*Kg[i]))/real(sum(conj(g_old[i]).*Kg_old[i]))
    end
    return
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

    Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()

    psiks = rand_BlochWavefunc( Ham )
    Etot_old = obj_function!( Ham, psiks, skip_ortho=true )

    g = zeros_BlochWavefunc( Ham )
    gt = zeros_BlochWavefunc( Ham )    
    Kg = zeros_BlochWavefunc( Ham )
    g_old = zeros_BlochWavefunc( Ham )
    Kg_old = zeros_BlochWavefunc( Ham )
    d = zeros_BlochWavefunc( Ham )
    d_old = zeros_BlochWavefunc( Ham )
    psic = zeros_BlochWavefunc( Ham )

    β = zeros(length(g))
    α = zeros(length(g))

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )


    α_t = 1e-5
    etot_conv_thr = 1e-7

    Nconverges = 0

    for iter = 1:50
        
        grad_obj_function!( Ham, psiks, g )
        precond_grad!( Ham, g, Kg )
        if iter > 1
            calc_beta_CG!( g, g_old, Kg, Kg_old, β )
        end

        d = -Kg + β .* d_old

        psic = psiks + α_t*d  # trial wavefunc

        # line minimization
        grad_obj_function!( Ham, psic, gt )
        calc_alpha_CG!( α_t, g, gt, d, α )

        psiks = psiks + α .* d

        Etot = obj_function!( Ham, psiks )
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

    end

end

main_CG()