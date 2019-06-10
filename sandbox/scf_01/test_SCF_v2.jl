using Printf
using Random
using LinearAlgebra

using PWDFT

include("alt2_KS_solve_SCF_spinpol.jl")

function create_Hamiltonian_Ni_fcc()
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Ni  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(6.65914911201) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Ni-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="LDA",
                       Nspin=2, meshk=[3,3,3], extra_states=4 )
    return Ham
end


function test_main()

    Ham = create_Hamiltonian_Ni_fcc()

    @time alt2_KS_solve_SCF!(
        Ham, etot_conv_thr=1e-6, NiterMax=100, betamix=0.5, update_psi="LOBPCG",
        mix_method="rpulay", check_rhoe=false,
        use_smearing=true, kT=1e-3
    )

end

test_main()
