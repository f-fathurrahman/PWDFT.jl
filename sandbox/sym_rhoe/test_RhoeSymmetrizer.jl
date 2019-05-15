using Printf
using LinearAlgebra
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")


function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Ns_=(32,32,32) )
end


function init_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.6839444516))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function init_Ham_atom_H()
    atoms = Atoms( xyz_string="""
            1

            H   0.1   0.0   0.0
            """, LatVecs=gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2_v2()
    atoms = Atoms( xyz_string=
        """
        2
            
        H  0.0  0.0  0.0
        H  1.5  0.0  0.0
        """, in_bohr=true, LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function main()

    Ham = init_Ham_Si_fcc()
    #Ham = init_Ham_GaAs()
    #Ham = init_Ham_atom_H()
    #Ham = init_Ham_H2()
    #Ham = init_Ham_H2_v2()

    pw = Ham.pw
    atoms = Ham.atoms
    LatVecs = pw.LatVecs

    sym_info = Ham.sym_info
    println(sym_info)

    rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    println("Ngs = ", rhoe_symmetrizer.Ngs)

    Random.seed!(1234)

    dVol = pw.CellVolume/prod(pw.Ns)
    psiks = rand_BlochWavefunc(Ham)
    Rhoe = calc_rhoe(Ham, psiks)
    println("integ Rhoe = ", sum(Rhoe)*dVol)

    Rhoe_sym = copy(Rhoe)
    symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe_sym )
    println("integ Rhoe_sym = ", sum(Rhoe_sym)*dVol)

    # Make sure that they are different
    println("diff sum abs = ", sum(abs.(Rhoe_sym - Rhoe)))

end

main()
