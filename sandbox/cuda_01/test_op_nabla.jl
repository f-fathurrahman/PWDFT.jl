using Test

include("PWDFT_cuda.jl")

function main()
    atoms = Atoms(xyz_string=
    """
    1

    H   0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc )
    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    @views ρ = Rhoe[:,1]
    ∇ρ = op_nabla( Ham.pw, ρ )
    @time CuArrays.@sync ∇ρ = op_nabla( Ham.pw, ρ )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )
    Rhoe_ = collect(Rhoe)

    @views ρ_ = Rhoe_[:,1]

    ∇ρ_cpu = op_nabla( Ham_.pw, Rhoe_[:,1] )
    @time ∇ρ_cpu = op_nabla( Ham_.pw, Rhoe_[:,1] )

    ∇ρ_gpu = collect(∇ρ)

    idir = 2
    for i in 1:5
        @printf("%18.10f %18.10f\n", ∇ρ_gpu[idir,i], ∇ρ_cpu[idir,i])
    end

    @test ∇ρ_gpu ≈ ∇ρ_cpu

    println("integ rhoe = ", sum(Rhoe)*Ham.pw.CellVolume/prod(Ham.pw.Ns))

    println("Several types: ")
    println( typeof(Rhoe_) )
    println( size(Rhoe_) )
    println( typeof(ρ) )
    println( typeof(ρ_) )
    println()

    println("Pass here")
end

main()
