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
    Si   0.0  1.51  0.0
    Si   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(5.0))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 5.0

    Nspin = 1

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    E_kin = calc_E_kin( Ham, psiks )

    E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham )

    E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )

    energies = calc_energies( Ham, psiks )

    @printf("E_kin     = %18.10f\n", E_kin)
    @printf("E_Ps_loc  = %18.10f\n", E_Ps_loc)
    @printf("E_Hartree = %18.10f\n", E_Hartree)
    @printf("E_xc      = %18.10f\n", E_xc)
    @printf("E_Ps_nloc = %18.10f\n", E_Ps_nloc)

    println(energies)

    #
    # Compare with CPU calculation
    #
    Nkspin = length(psiks)
    Ham_cpu = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, Nspin=Nspin )
    psiks_cpu = BlochWavefunc(undef, Nkspin)
    for i in 1:Nkspin
        psiks_cpu[i] = collect(psiks[i])
    end

    Rhoe_cpu = calc_rhoe( Ham_cpu, psiks_cpu )
    update!( Ham_cpu, Rhoe_cpu )

    E_kin_cpu = calc_E_kin( Ham_cpu, psiks_cpu )

    E_Ps_loc_cpu, E_Hartree_cpu, E_xc_cpu = calc_E_local( Ham_cpu )

    E_Ps_nloc_cpu = calc_E_Ps_nloc( Ham_cpu, psiks_cpu )

    energies_cpu = calc_energies( Ham_cpu, psiks_cpu )

    @printf("E_kin_cpu     = %18.10f\n", E_kin_cpu)
    @printf("E_Ps_loc_cpu  = %18.10f\n", E_Ps_loc_cpu)
    @printf("E_Hartree_cpu = %18.10f\n", E_Hartree_cpu)
    @printf("E_xc_cpu      = %18.10f\n", E_xc_cpu)
    @printf("E_Ps_nloc_cpu = %18.10f\n", E_Ps_nloc_cpu)

    @test E_kin ≈ E_kin_cpu
    @test E_Ps_loc ≈ E_Ps_loc_cpu
    @test E_Hartree ≈ E_Hartree_cpu
    @test E_xc ≈ E_xc_cpu
    @test E_Ps_nloc ≈ E_Ps_nloc_cpu

    println("Pass here")
end

main()