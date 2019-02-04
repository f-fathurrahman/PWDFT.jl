using LinearAlgebra
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_E_NN_v2.jl")

#=
ewald = 3.58126919 Ry = 1.790634595 Ha
=#
function test_N2()
    atoms = Atoms(xyz_file=joinpath(DIR_STRUCTURES,"N2.xyz"))

    LatVecs = gen_lattice_sc(16.0)
    Zvals = [5.0] # XXX hardwired

    E_NN = calc_E_NN_v2( LatVecs, atoms, Zvals )

    E_NN_ref = 1.790634595
    diffE = abs(E_NN - E_NN_ref)

    @printf("N2    : E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
ewald contribution = 0.62633998 Ry = 0.31316999 Ha
"""
function test_H2()
    atoms = Atoms(xyz_file=joinpath(DIR_STRUCTURES,"H2.xyz"))

    LatVecs = gen_lattice_sc(16.0)
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN_v2( LatVecs, atoms, Zvals )

    E_NN_ref = 0.31316999 
    diffE = abs(E_NN - E_NN_ref)

    @printf("H2    : E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
PWSCF result = -0.04393722 Ry = -0.02196861 Ha
"""
function test_LiH()
    atoms = Atoms(xyz_file=joinpath(DIR_STRUCTURES,"LiH.xyz"))

    LatVecs = gen_lattice_sc(16.0)
    Zvals = [1.0,1.0] # XXX hardwired

    E_NN = calc_E_NN_v2( LatVecs, atoms, Zvals )
    
    E_NN_ref = -0.02196861
    diffE = abs(E_NN - E_NN_ref)

    @printf("LiH   : E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end

"""
PWSCF result = -0.17733109 Ry = -0.088665545 Ha
"""
function test_H()
    atoms = Atoms(xyz_file=joinpath(DIR_STRUCTURES,"H.xyz"))

    LatVecs = gen_lattice_sc(16.0)
    Zvals = [1.0] # XXX hardwired

    E_NN = calc_E_NN_v2( LatVecs, atoms, Zvals )

    E_NN_ref = -0.088665545
    diffE = abs(E_NN - E_NN_ref)

    @printf("H     : E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)

end


function test_Si_fcc()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    atoms.Zvals = [4.0]  # hardwired

    E_NN = calc_E_NN_v2( atoms )

    E_NN_ref = -16.79585054*0.5

    diffE = abs(E_NN - E_NN_ref)

    @printf("Si_fcc: E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)


end


function test_GaAs()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

    atoms.Zvals = [3.0, 5.0]  # hardwired

    E_NN = calc_E_NN_v2( atoms )

    E_NN_ref = -16.84241132*0.5

    diffE = abs(E_NN - E_NN_ref)

    @printf("GaAs  : E_NN, ref, diff: %18.10f %18.10f %18.10e\n", E_NN, E_NN_ref, diffE)


end

test_H()
test_H2()
test_LiH()
test_N2()
test_Si_fcc()
test_GaAs()