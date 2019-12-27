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

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks = rand_CuBlochWavefunc( Ham )
    
    prec_psi = Kprec( 1, Ham.pw, psiks[1] )

    # Compare with CPU

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    Nkspin = length(psiks)
    psiks_ = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_[i] = collect(psiks[i])
    end
    prec_psi_cpu = Kprec( 1, Ham_.pw, psiks_[1] )

    prec_psi_gpu = collect( prec_psi )

    ist = 1
    println("Some prec_psi:")
    for ig in 1:5
        v1 = prec_psi_cpu[ig,ist]
        v2 = prec_psi_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    @test prec_psi_cpu ≈ prec_psi_gpu

    println("Pass here")
end


main()
