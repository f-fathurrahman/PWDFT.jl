include("PWDFT_cuda.jl")

function main()
    atoms = Atoms(xyz_string=
    """
    3

    H   0.0  0.0  0.0
    H   0.0  1.5  0.0
    H   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc )

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)
    println("integ Rhoe = ", sum(Rhoe)*dVol)

    println("Pass here")
end

main()