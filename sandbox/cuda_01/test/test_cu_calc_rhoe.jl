using Test
using Random
using Printf

using PWDFT
using PWDFT_cuda

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string=
    """
    3

    H   0.0  0.0   0.0
    H   0.0  1.51   0.0
    H   1.5  0.0   0.0
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 5.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )
    
    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks_cpu = BlochWavefunc(undef, Nkspin)
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
    for i in 1:10
        @printf("%8d %18.10f %18.10f\n", i, Rhoe_cpu[i,1], Rhoe_gpu[i,1])
    end

    println(typeof(Rhoe))
    println(typeof(psiks))

    println(typeof(Rhoe_cpu))
    println(typeof(Rhoe_gpu))

    # Generally this will not work because ortho_sqrt and ortho_gram_schmidt will not
    # give the same result
    @test Rhoe_cpu ≈ Rhoe_gpu

    println("Pass here")
end

main()