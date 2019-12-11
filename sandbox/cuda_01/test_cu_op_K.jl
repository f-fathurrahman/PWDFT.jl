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
    Kpsi = op_K( Ham, psiks[1] )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )
    Nkspin = length(psiks)
    psiks_ = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_[i] = collect(psiks[i])
    end
    Kpsi_cpu = op_K( Ham_, psiks_[1] )

    Kpsi_gpu = collect( Kpsi )

    println("diff = ", sum(Kpsi_cpu - Kpsi_gpu) )

    println("Pass here")
end


main()
