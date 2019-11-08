using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("grad_obj_function.jl")
include("create_Ham.jl")

function create_Ham_Si_fcc( ; xcfunc="VWN", Nspin=1 )

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]

    ecutwfc = 15.0
    if xcfunc == "PBE"
        if Nspin == 2
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], xcfunc="PBE" )
        else
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], xcfunc="PBE", Nspin=2, extra_states=4 )
        end
    else
        if Nspin == 2
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2, extra_states=4 )
        else
            return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
        end
    end
end

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
        if denum != 0.0
            #α[i] = abs( α_t*real(sum(conj(g[i]).*d[i]))/denum )
            α[i] = α_t*real(sum(conj(g[i]).*d[i]))/denum
        else
            α[i] = 0.0
        end
    end
    return
end


include("quadratic_line_minimizer.jl")

function main_CG()

    Random.seed!(1234)

    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    Ham = create_Ham_Si_fcc()  # pathological for current implementation of quadratic_line_minimizer

    psiks = rand_BlochWavefunc( Ham )
    Etot_old = obj_function!( Ham, psiks, skip_ortho=true )

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64,Npoints,Nspin)
    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )
    #
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, psiks, verbose=false, NiterMax=20 )


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

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    α_t = 1.0
    etot_conv_thr = 1e-6

    Nconverges = 0

    for iter = 1:20

        grad_obj_function!( Ham, psiks, g )
        precond_grad!( Ham, g, Kg )
        if iter > 1
            calc_beta_CG!( g, g_old, Kg, Kg_old, β )
        end
        println("β = ", β)
        d = -Kg + β .* d_old

        # line minimization
        #psic = psiks + α_t*d  # trial wavefunc
        #grad_obj_function!( Ham, psic, gt )
        #calc_alpha_CG!( α_t, g, gt, d, α )

        α_t = quadratic_line_minimizer!(Ham, psiks, d, Etot_old, g, α_t, α)

        println("α = ", α)

        psiks = psiks + α .* d

        Etot = obj_function!( Ham, psiks )
        diffE = Etot_old - Etot
        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        if abs(diffE) < etot_conv_thr
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

@time main_CG()
#@time main_CG()
