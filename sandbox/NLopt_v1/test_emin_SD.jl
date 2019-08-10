using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("obj_function.jl")
include("grad_obj_function.jl")
include("create_Ham.jl")

# steepest descent
function main_SD(; etot_conv_thr=1e-6)

    Random.seed!(1234)

    #Ham = create_Ham_H2()
    Ham = create_Ham_H_atom()

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
    _ = diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )

    g = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )

    α_t = 3e-5

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    Nconverges = 0

    for iter = 1:50

        grad_obj_function!( Ham, psiks, g )
        precond_grad!( Ham, g, Kg )
        
        psiks = psiks - α_t*Kg
        
        Etot = obj_function!( Ham, psiks )
        
        diffE = Etot_old - Etot
        @printf("%8d %18.10f %18.10e\n", iter, Etot, diffE)

        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nEmin_SD is converged in iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end


        Etot_old = Etot
    end

end

main_SD()