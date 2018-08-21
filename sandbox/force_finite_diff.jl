using LinearAlgebra
using Printf
using PWDFT

function calc_force_finite_diff( atoms::Atoms, pspfiles, ecutwfc )

    pos_orig = copy( atoms.positions )
    
    Natoms = atoms.Natoms
    force = zeros(3,Natoms)

    Δ = 0.005
    for ia = 1:Natoms
    for ii = 1:3
        atoms.positions[ii,ia] = pos_orig[ii,ia] + Δ
        Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
        #
        KS_solve_SCF!( Ham, mix_method="anderson" )
        Etot1 = sum(Ham.energies)
        println("Etot1 = ", Etot1)
        #
        atoms.positions[ii,ia] = pos_orig[ii,ia] - Δ
        Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
        #
        KS_solve_SCF!( Ham, mix_method="anderson" )
        Etot2 = sum(Ham.energies)
        println("Etot2 = ", Etot2)
        #
        force[ii,ia] = (Etot1 - Etot2) / (2Δ)
    end
    end

    return force

end

function test_main()

    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """
    )
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)

    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0

    force = calc_force_finite_diff(atoms, pspfiles, ecutwfc)

    println("force = ")
    println(force)

end

test_main()
