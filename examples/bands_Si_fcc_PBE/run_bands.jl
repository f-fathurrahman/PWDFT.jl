using Printf
using LinearAlgebra
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pbe_gth")

include("dump_bandstructure.jl")
include("gen_kpath.jl")

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(10.2631))
    
    pspfiles = [joinpath(DIR_PSP,"Si-q4.gth")]
    ecutwfc = 15.0

    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4, xcfunc="PBE" )

    KS_solve_SCF!( Ham, mix_method="rpulay" )

    # Band structure calculation
    kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "WLGXWK", "fcc", Δk=0.05 )

    # New pw
    pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)
    Ham.pw = pw

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=1, Nkpt=kpoints.Nkpt,
                               Nstates_empty=1 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( atoms, pw, Ham.pspots )

    psiks = rand_BlochWavefunc( Ham )
    evals = zeros(Float64,Nstates,Nkspin)

    k = Ham.pw.gvecw.kpoints.k
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        #
        @printf("\nispin = %d, ik = %d, ikspin=%d, Ngw = %d\n", ispin, ik, ikspin, Ngw[ik])
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        evals[:,ikspin], psiks[ikspin] =
        diag_LOBPCG( Ham, psiks[ikspin], verbose_last=true )
    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_bands.dat" )

end

main()
