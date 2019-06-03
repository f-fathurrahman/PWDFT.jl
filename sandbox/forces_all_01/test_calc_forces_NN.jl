using LinearAlgebra
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_forces_NN.jl")

function test_H2()

    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs=gen_lattice_sc(16.0))

    atoms.Zvals = [1.0]  # hardwired

    println("")
    F_NN = calc_forces_NN( atoms )*2.0

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_NN = calc_forces_NN_finite_diff( atoms )*2.0

    println("")
    println("Using finite difference")
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

end



function test_GaAs()

    LatVecs = zeros(3,3)
    LatVecs[:,1] = [0.5, 0.5, 0.0]
    LatVecs[:,2] = [0.5, 0.0, 0.5]
    LatVecs[:,3] = [0.0, 0.5, 0.5]
    LatVecs = LatVecs*5.6537*ANG2BOHR

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

#    atoms = Atoms(xyz_string_frac=
#        """
#        2
#
#        Ga  0.0  0.0  0.0
#        As  0.25  0.35  0.25
#        """, in_bohr=true, LatVecs=LatVecs)

    atoms.Zvals = [3.0, 5.0]  # hardwired

    println("")
    F_NN = calc_forces_NN( atoms )*2.0

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_NN = calc_forces_NN_finite_diff( atoms )*2.0

    println("")
    println("Using finite difference")
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

end


function test_Si_fcc()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.15  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    atoms.Zvals = [4.0]  # hardwired

    println("")
    F_NN = calc_forces_NN( atoms )*2.0  # convert to Ry

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_NN = calc_forces_NN_finite_diff( atoms )*2.0

    println("")
    println("Using finite difference")
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

end

test_H2()
test_GaAs()
test_Si_fcc()

