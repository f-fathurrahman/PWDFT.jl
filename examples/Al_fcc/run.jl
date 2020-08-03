using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

import InteractiveUtils
InteractiveUtils.versioninfo()
import Dates
println("Now = ", Dates.now())

function main()

    Random.seed!(1234)

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(7.6525970200) )
    
    write_xsf("TEMP_Al.xsf", atoms)

    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )
    println(Ham)

    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2, use_smearing=true )
    
end

@time main()
@time main()
