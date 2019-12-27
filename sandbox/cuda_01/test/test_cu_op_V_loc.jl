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

    H   0.0  0.0  0.0
    H   0.0  1.51  0.0
    H   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    Vpsi = op_V_loc( Ham, psiks[1] )

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

    update!( Ham_cpu, Rhoe_cpu )

    Vpsi_cpu = op_V_loc( Ham_cpu, psiks_cpu[1] )

    #
    # Compare
    #

    Vpsi_gpu = collect(Vpsi)

    ist = 1
    println("Some Vpsi")
    for ig in 1:5
        v1 = Vpsi_cpu[ig,ist]
        v2 = Vpsi_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    @test Vpsi_cpu ≈ Vpsi_gpu

    println("Pass here")
end

main()