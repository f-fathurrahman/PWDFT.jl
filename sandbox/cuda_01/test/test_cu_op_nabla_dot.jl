using Test
using Random
using Printf

using CuArrays
using PWDFT
using PWDFT_cuda

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()
    atoms = Atoms(xyz_string=
    """
    1

    H   0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    @views ρ = Rhoe[:,1]
    ∇ρ = op_nabla( Ham.pw, ρ )
    ∇2ρ = op_nabla_dot( Ham.pw, ∇ρ )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    Rhoe_ = collect(Rhoe)

    ∇ρ_cpu = op_nabla( Ham_.pw, Rhoe_[:,1] )
    ∇2ρ_cpu = op_nabla_dot( Ham_.pw, ∇ρ_cpu )

    ∇2ρ_gpu = collect(∇2ρ)

    for i in 1:5
        @printf("%18.10f %18.10f\n", ∇2ρ_gpu[i], ∇2ρ_cpu[i])
    end

    @test ∇2ρ_gpu ≈ ∇2ρ_cpu

    println("integ rhoe = ", sum(Rhoe)*Ham.pw.CellVolume/prod(Ham.pw.Ns))

    println("Pass here")
end

main()
