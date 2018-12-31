using Printf
using PWDFT

include("ABINIT.jl")

function test_CuSO4()
    # initialize atoms and Hamiltonian
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )

    write_abinit( Ham, prefix_psp="../_compare/abinit_psp/",
                       prefix="../TEMP_ABINIT/" )
end
#test_CuSO4()


function test_H2()
    # initialize atoms and Hamiltonian
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )

    write_abinit( Ham, prefix_psp="../_compare/abinit_psp/",
                       prefix="../TEMP_ABINIT/" )
end
#test_H2()

function test_Al_fcc()

    atoms = Atoms( xyz_string_frac=
    """
    1

    Al  0.0  0.0  0.0
    """, in_bohr=true,
    LatVecs = gen_lattice_fcc(7.6525970200) )
    pspfiles = ["../pseudopotentials/pade_gth/Al-q3.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )

    write_abinit( Ham, use_smearing=true,
                  prefix_psp="../_compare/abinit_psp/",
                  prefix="../TEMP_ABINIT/" )
end
#test_Al_fcc()

function test_read_abinit_etotal()
    read_abinit_etotal("../TEMP_ABINIT/LOG1")
end
#test_read_abinit_etotal()

function test_Pt_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE",
                        meshk=[3,3,3], extra_states=4 )
    
    write_abinit( Ham, use_smearing=true,
                  prefix_psp="../_compare/abinit_psp/",
                  prefix="../TEMP_ABINIT/" )
end
test_Pt_fcc()