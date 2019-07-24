using Printf
using LinearAlgebra
using Random
using PWDFT
import Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include("../common/dump_bandstructure.jl")
include("../common/gen_kpath.jl")

function main()

    Random.seed!(1234)

    Ham = Serialization.deserialize("Ham.data")

    atoms = Ham.atoms
    ecutwfc = Ham.pw.ecutwfc
    Nspin = Ham.electrons.Nspin

    # Band structure calculation
    kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "W-L-G-X-W-K", "fcc", Δk=0.05 )
    #kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "L-G-X-G1", "fcc", Δk=0.05 )
    #kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "G-X-W-K-G-L-U-W-L-K", "fcc", Δk=0.05 )

    # New pw
    pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)
    Ham.pw = pw

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                               Nstates_empty=4 )

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
