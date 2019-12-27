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

    Si   0.0  0.0  0.0
    Si   0.0  1.51  0.0
    Si   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 5.0

    Nspin = 1

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )

    g = calc_grad( Ham, psiks[1] )

    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    psiks_cpu = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_cpu[i] = collect(psiks[i])
    end

    Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )
    update!( Ham_cpu, Rhoe_cpu )

    g_cpu = calc_grad( Ham_cpu, psiks_cpu[1] )

    g_gpu = collect(g)

    ist = 1
    println("Some g:")
    for ig in 1:5
        v1 = g_cpu[ig,ist]
        v2 = g_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    @test g_cpu ≈ g_gpu

    println("Pass here")
end

main()