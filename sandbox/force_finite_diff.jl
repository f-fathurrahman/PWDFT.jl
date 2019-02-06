using LinearAlgebra
using Printf
using PWDFT

function calc_force_finite_diff( atoms::Atoms, pspfiles, ecutwfc )

    pos_orig = copy(atoms.positions)
    
    Natoms = atoms.Natoms
    forces = zeros(3,Natoms)

    forces_NN = zeros(3,Natoms)
    forces_Ha = zeros(3,Natoms)
    forces_Ps_loc = zeros(3,Natoms)
    forces_XC = zeros(3,Natoms)
    forces_Kinetic = zeros(3,Natoms)
    forces_Ps_nloc = zeros(3,Natoms)

    Δ = 0.005
    for ia = 1:Natoms
    for ii = 1:3
        
        atoms.positions[:,:] = copy(pos_orig)

        atoms.positions[ii,ia] = pos_orig[ii,ia] + 0.5*Δ

        println("")
        println("positions in angstrom")
        for ia = 1:Natoms
            @printf("%s %18.10f %18.10f %18.10f\n", atoms.atsymbs[ia],
                atoms.positions[1,ia]/ANG2BOHR,
                atoms.positions[2,ia]/ANG2BOHR,
                atoms.positions[3,ia]/ANG2BOHR)
        end
        
        Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
        #
        KS_solve_SCF!( Ham, mix_method="rpulay", ETOT_CONV_THR=1e-8 )
        Eplus = sum(Ham.energies)
        #
        atoms.positions[ii,ia] = pos_orig[ii,ia] - 0.5*Δ
        
        println("")
        println("positions in angstrom")
        for ia = 1:Natoms
            @printf("%s %18.10f %18.10f %18.10f\n", atoms.atsymbs[ia],
                atoms.positions[1,ia]/ANG2BOHR,
                atoms.positions[2,ia]/ANG2BOHR,
                atoms.positions[3,ia]/ANG2BOHR)
        end

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
        #
        KS_solve_SCF!( Ham, mix_method="rpulay", ETOT_CONV_THR=1e-8 )
        Eminus = sum(Ham.energies)
        #
        forces[ii,ia] = -(Eplus - Eminus)/Δ

        println("")
        @printf("ia = %d, idir = %d\n", ia, ii)
        @printf("Eplus  = %18.10f\n", Eplus)
        @printf("Eminus = %18.10f\n", Eminus)
        @printf("diff   = %18.10f\n", Eplus - Eminus)
        @printf("F      = %18.10f\n", forces[ii,ia])
    end
    end

    return forces

end

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

function test_main()

    # Atoms
    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs = gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    forces = calc_force_finite_diff(atoms, pspfiles, ecutwfc)*2.0

    println("")
    println("forces (in Ry/au) = ")
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                forces[1,ia], forces[2,ia], forces[3,ia] )
    end

end

test_main()
