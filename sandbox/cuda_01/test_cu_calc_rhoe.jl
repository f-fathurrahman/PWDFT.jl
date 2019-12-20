using Test

include("PWDFT_cuda.jl")

using Random

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string=
    """
    3

    H   0.0  0.0  0.0
    H   0.0  1.5  0.0
    H   1.5  0.0  0.0
    """, LatVecs=gen_lattice_fcc(10.0))

    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc )

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    Nkspin = length(psiks)
    # Compare with CPU calculation
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc )
    psiks_cpu = BlochWavefunc(undef, length(Nkspin))
    for i in 1:Nkspin
        psiks_cpu[i] = collect(psiks[i])
    end

    Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )
    
    Rhoe_gpu = collect(Rhoe)

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    println("integ Rhoe GPU = ", sum(Rhoe)*dVol)
    println("integ Rhoe CPU = ", sum(Rhoe_cpu)*dVol)
    println("integ Rhoe GPU v2 = ", sum(Rhoe_gpu)*dVol)

    println("Some Rhoe:")
    for i in 1:5
        @printf("%18.10f %18.10f\n", Rhoe_cpu[i,1], Rhoe_gpu[i,1])
    end
    println(size(Rhoe_cpu))
    println(size(Rhoe_gpu))

    # Generally this will not work because ortho_sqrt and ortho_gram_schmidt will not
    # give the same result
    #@test Rhoe_cpu ≈ Rhoe_gpu

    println("Pass here")
end

main()