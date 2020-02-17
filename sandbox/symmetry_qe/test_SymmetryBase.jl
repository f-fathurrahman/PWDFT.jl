using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("SymmetryBase.jl")

function create_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.25
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    # FIXME: Zvals is not set
    return atoms
end

function main()
    atoms = create_Si_fcc()

    nrot, s, ft, sname = find_symm_bravais_latt( atoms.LatVecs )

    sym = zeros(Bool, 48)
    irt = zeros(Int64, 48, atoms.Natoms)

    sgam_at!( atoms, sym, false, nrot, s, ft, sname, irt )

    Nsyms = copy_sym!( nrot, sym, s, ft, sname, irt )
    println("Nsyms = ", Nsyms)

    for i in 1:Nsyms
        println()
        println("symmetry operation: ", sname[i])
        display(s[:,:,i]); println()
        println("fractional translations: ")
        println(ft[:,i])
    end
    println("irt = ", irt[1:Nsyms,:])
end

main()