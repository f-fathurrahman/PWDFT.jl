using Printf
using Random
using LinearAlgebra

using PWDFT

function create_Hamiltonian_H_atom()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """)
    atoms.LatVecs = gen_lattice_sc(16.0)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    return Ham
end

function create_Hamiltonian_H2()
    # Atoms
    atoms = init_atoms_xyz("../structures/H2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    return Ham
end

function create_Hamiltonian_N2()
    # Atoms
    atoms = init_atoms_xyz("../structures/N2.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    return Ham
end


function create_Hamiltonian_O2()
    # Atoms
    atoms = init_atoms_xyz("../structures/O2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/O-q6.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, extra_states=1 )
    Ham.electrons.Focc[:,1] = [2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0]

    println(Ham)

    return Ham
end

function create_Hamiltonian_O2_spinpol()
    # Atoms
    atoms = init_atoms_xyz("../structures/O2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/O-q6.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc_Ry*0.5, extra_states=4, Nspin=2
    )
    Ham.electrons.Focc[:,1] = [1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 0.0, 0.0, 0.0]
    Ham.electrons.Focc[:,2] = [1.0, 1.0, 1.0, 1.0, 1.0,
                               0.0, 0.0, 0.0, 0.0, 0.0]                               

    println(Ham)

    return Ham
end


function create_Hamiltonian_N2()
    # Atoms
    atoms = init_atoms_xyz("../structures/N2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, extra_states=1 )

    return Ham
end

function create_Hamiltonian_Co_atom()
    atoms = Atoms(xyz_string="""
    1

    Co   0.0   0.0   0.0
    """, LatVecs=gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/Co-q9.gth"]
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc_Ry*0.5, Nspin=2, extra_states=3
        )

    # set manually
    #Ham.electrons.Focc[:,1] = [1.0, 1.0, 1.0, 1.0, 0.5, 0.0]
    #Ham.electrons.Focc[:,2] = [1.0, 1.0, 1.0, 1.0, 0.5, 0.0]
    #println("\nFocc is set manually\n")
    #println(Ham.electrons)

    return Ham
end

function create_Hamiltonian_Ni_atom()
    atoms = Atoms(xyz_string="""
    1

    Ni   0.0   0.0   0.0
    """, LatVecs=gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/Ni-q10.gth"]
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc_Ry*0.5, Nspin=2, extra_states=3
        )

    # set manually
    Ham.electrons.Focc[:,1] = [1.0, 1.0, 1.0, 1.0, 1.0,
                               0.0, 0.0, 0.0]
    Ham.electrons.Focc[:,2] = [1.0, 1.0, 1.0, 1.0, 0.5,
                               0.5, 0.0, 0.0]

    println(Ham)

    return Ham
end


function test_main()

    #Ham = create_Hamiltonian_H_atom()
    Ham = create_Hamiltonian_H2()
    #Ham = create_Hamiltonian_N2()
    #Ham = create_Hamiltonian_O2()
    #Ham = create_Hamiltonian_Co_atom()
    #Ham = create_Hamiltonian_Ni_atom()
    #Ham = create_Hamiltonian_O2_spinpol()

    # Solve the KS problem
    @time KS_solve_SCF!(
        Ham, etot_conv_thr=1e-6, NiterMax=100, betamix=0.7, update_psi="LOBPCG",
        mix_method="pulay"
    )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    if Ham.electrons.Nspin == 2
        for ist = 1:Nstates
            @printf("%8d : %18.10f = %18.10f eV -- %18.10f = %18.10f eV\n", ist, ebands[ist,1], ebands[ist,1]*Ry2eV*2, ebands[ist,2], ebands[ist,2]*Ry2eV*2)
        end
    else
        for ist = 1:Nstates
            @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist,1], ebands[ist,1]*Ry2eV*2)
        end
    end
    
    println("\nTotal energy components")
    println(Ham.energies)
end

test_main()
