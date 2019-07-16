using PWDFT
using Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

function test_pwgrid()

    pw1 = PWGrid( 15.0, gen_lattice_fcc(10.0) )
    println(pw1)
    serialize( "pw.data", pw1 )

    pw2 = deserialize("pw.data")
    println(pw2)

end
#test_pwgrid()


function create_Ham_Al_fcc()
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(7.6525970200) )
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )
end

function test_hamiltonian()
    Ham1 = create_Ham_Al_fcc()
    serialize("Ham.data", Ham1)

    Ham1 = 0

    Ham2 = deserialize("Ham.data")
    println(Ham2)
end
test_hamiltonian()
