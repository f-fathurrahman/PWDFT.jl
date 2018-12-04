using Printf
using PWDFT

function my_scf!( Ham::Hamiltonian; NiterMax=100, betamix=0.2 )

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nspin*Nkpt
    Nstates = Ham.electrons.Nstates
    atoms = Ham.atoms
    pspots = Ham.pspots
    CellVolume = Ham.pw.CellVolume

    psiks = rand_BlochWavefunc(Ham)
    
    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_new = zeros(Float64,Npoints,Nspin)

    Rhoe[:,:] = calc_rhoe(Ham, psiks)

    update!(Ham, Rhoe)
    Ham.energies.NN = calc_E_NN(atoms)
    Ham.energies.PspCore = calc_PspCore_ene(atoms, pspots, CellVolume)

    evals = zeros(Nstates,Nkspin)

    Etot_old = 0.0

    @printf("\n")
    @printf("SCF iteration starts\n")
    @printf("\n")

    for iterSCF = 1:NiterMax
        evals = diag_LOBPCG!( Ham, psiks )
        Rhoe_new[:,:] = calc_rhoe( Ham, psiks )
        Rhoe[:,:] = betamix*Rhoe_new + (1-betamix)*Rhoe
        update!(Ham, Rhoe)
        Ham.energies = calc_energies(Ham, psiks)
        Etot = sum(Ham.energies)
        diffEtot = abs(Etot - Etot_old)
        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        if diffEtot <= 1e-6
            @printf("SCF is converged in %d iterations\n", iterSCF)
            return
        end
        Etot_old = Etot
    end
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end

function test_my_scf()

    DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
    atoms = Atoms(
        xyz_file=joinpath(DIR_PWDFT, "structures/H2.xyz"),
        LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, 15.0 )
    println(Ham)

    my_scf!(Ham, betamix=0.5)

    println()
    println("Kohn-Sham energy components\n")
    println(Ham.energies)
end

test_my_scf()
