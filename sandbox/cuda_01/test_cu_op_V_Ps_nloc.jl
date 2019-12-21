using Test

include("PWDFT_cuda.jl")

function main()
    atoms = Atoms(xyz_string=
    """
    1

    H   0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks = rand_CuBlochWavefunc( Ham )

    println("NbetaNL = ", Ham.pspotNL.NbetaNL)
    @assert Ham.pspotNL.NbetaNL > 0

    Vpsi = op_V_Ps_nloc( Ham, psiks[1] )

    # Compare with CPU

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    Nkspin = length(psiks)
    psiks_ = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_[i] = collect(psiks[i])
    end
    Vpsi_cpu = op_V_Ps_nloc( Ham_, psiks_[1] )

    Vpsi_gpu = collect( Vpsi )

    ist = 1
    println("Some Vpsi:")
    for ig in 1:5
        v1 = Vpsi_cpu[ig,ist]
        v2 = Vpsi_gpu[ig,ist]
        @printf("%18.10f + im*%18.10f, %18.10f + im*%18.10f\n", real(v1), imag(v1), real(v2), imag(v2) )
    end

    #@test Kpsi_cpu ≈ Kpsi_gpu

    println("Pass here")
end


main()