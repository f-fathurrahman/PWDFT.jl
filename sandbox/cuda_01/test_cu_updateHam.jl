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
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 5.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc )
    psiks_cpu = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_cpu[i] = collect(psiks[i])
    end

    Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )

    update!( Ham_cpu, Rhoe_cpu )

    println("Pass here")
end

main()