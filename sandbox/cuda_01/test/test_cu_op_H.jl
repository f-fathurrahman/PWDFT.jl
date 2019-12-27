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
    H    0.0  1.51  0.0
    Si   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth"), 
                joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    Hpsi = op_H( Ham, psiks[1] )

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

    Hpsi_cpu = op_H( Ham_cpu, psiks_cpu[1] )

    #
    # Compare
    #

    Hpsi_gpu = collect(Hpsi)

    ist = 1
    println("Some Hpsi")
    for ig in 1:5
        v1 = Hpsi_cpu[ig,ist]
        v2 = Hpsi_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    @test Hpsi_cpu ≈ Hpsi_gpu

    println("Pass here")
end

main()