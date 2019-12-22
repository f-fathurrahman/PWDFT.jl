using Test

include("PWDFT_cuda.jl")

using Random

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string=
    """
    3

    Si   0.0  0.0  0.0
    Si   0.0  1.51  0.0
    Si   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc = 5.0

    Nspin = 1

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )

    #
    # Compare with CPU calculation
    #
    #Nkspin = length(psiks)
    #Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

    #println("Nsyms = ", Ham_cpu.sym_info.Nsyms)

    #psiks_cpu = BlochWavefunc(undef, Nkspin)
    #for i in 1:Nkspin
    #    psiks_cpu[i] = collect(psiks[i])
    #end

    #Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )

    #update!( Ham_cpu, Rhoe_cpu )

    println("Pass here")
end

main()