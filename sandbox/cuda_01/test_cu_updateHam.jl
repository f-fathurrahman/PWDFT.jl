using Test

include("PWDFT_cuda.jl")

using Random

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string=
    """
    3

    H   0.0  0.0  0.0
    H   0.0  1.51  0.0
    H   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 5.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    println("Nsyms = ", Ham_cpu.sym_info.Nsyms)

    psiks_cpu = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_cpu[i] = collect(psiks[i])
    end

    Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )

    update!( Ham_cpu, Rhoe_cpu )

    #
    # Compare
    #

    ispin = 1

    Rhoe_gpu = collect(Rhoe)
    println("Some Rhoe")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", Rhoe_cpu[ip,1], Rhoe_gpu[ip,1])
    end

    VH_cpu = Ham_cpu.potentials.Hartree
    VH_gpu = collect(Ham.potentials.Hartree)
    println("Some V Hartree")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", VH_cpu[ip,1], VH_gpu[ip,1])
    end

    Vps_loc_cpu = Ham_cpu.potentials.Ps_loc
    Vps_loc_gpu = collect(Ham.potentials.Ps_loc)
    println("Some V Ps_loc")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", Vps_loc_cpu[ip,1], Vps_loc_gpu[ip,1])
    end

    Vxc_cpu = Ham_cpu.potentials.XC
    Vxc_gpu = collect(Ham.potentials.XC)

    println("Some Vxc")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", Vxc_cpu[ip,1], Vxc_gpu[ip,1])
    end

    println("Pass here")
end

main()