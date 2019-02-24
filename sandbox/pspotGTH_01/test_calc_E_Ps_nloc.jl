using Random
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function test_main()

    # Atoms
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"PtO.xyz"),
        LatVecs=gen_lattice_sc(20.0)
    )

    pspfiles = [joinpath(DIR_PSP,"Pt-q10.gth"),
                joinpath(DIR_PSP,"O-q6.gth")]

    ecutwfc = 60.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=2 )
    
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)

    if Ham.pspotNL.NbetaNL > 0
        E_ps_NL = calc_E_Ps_nloc( Ham, psiks )
        println("E ps NL = ", E_ps_NL)
    end
    

end

test_main()
