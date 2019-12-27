using Test
using Random
using Printf

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
    Kpsi = op_K( Ham, psiks[1] )

    # Compare with CPU

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    Nkspin = length(psiks)
    psiks_ = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_[i] = collect(psiks[i])
    end
    Kpsi_cpu = op_K( Ham_, psiks_[1] )

    Kpsi_gpu = collect( Kpsi )


    ist = 1
    println("Some Kpsi:")
    for ig in 1:5
        v1 = Kpsi_cpu[ig,ist]
        v2 = Kpsi_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    @test Kpsi_cpu ≈ Kpsi_gpu

    println("Pass here")
end


main()
