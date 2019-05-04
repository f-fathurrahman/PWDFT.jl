using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
end


function do_calc()

    LATCONST = 10.2631
    Ham = init_Ham_Si_fcc( LATCONST, [3,3,3] )
    println(Ham)

    filname = "TEMP_vltot.xsf"
    write_xsf( filname, Ham.atoms )
    write_xsf_data3d_crystal( filname, Ham.atoms, Ham.pw.Ns, Ham.potentials.Ps_loc )

end

do_calc()