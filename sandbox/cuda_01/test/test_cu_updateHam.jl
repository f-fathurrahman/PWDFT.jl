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

    Nspin = 2

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

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
        @printf("%18.10f %18.10f\n", Rhoe_cpu[ip,ispin], Rhoe_gpu[ip,ispin])
    end

    @test Rhoe_cpu ≈ Rhoe_gpu

    VH_cpu = Ham_cpu.potentials.Hartree
    VH_gpu = collect(Ham.potentials.Hartree)
    println("Some V Hartree")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", VH_cpu[ip], VH_gpu[ip])
    end

    @test VH_cpu ≈ VH_gpu

    Vps_loc_cpu = Ham_cpu.potentials.Ps_loc
    Vps_loc_gpu = collect(Ham.potentials.Ps_loc)
    println("Some V Ps_loc")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", Vps_loc_cpu[ip], Vps_loc_gpu[ip])
    end

    @test Vps_loc_cpu ≈ Vps_loc_gpu

    Vxc_cpu = Ham_cpu.potentials.XC
    Vxc_gpu = collect(Ham.potentials.XC)

    println("Some Vxc")
    for ip = 1:5
        @printf("%18.10f %18.10f\n", Vxc_cpu[ip,ispin], Vxc_gpu[ip,ispin])
    end

    @test Vxc_cpu ≈ Vxc_gpu

    println("Pass here")
end

main()