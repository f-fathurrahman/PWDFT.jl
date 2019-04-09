using Printf
using PWDFT

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function calc_v_of_0( Ham::Hamiltonian )

    atoms = Ham.atoms
    pspots = Ham.pspots

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms

    v_of_0 = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        rloc = psp.rlocal
        zval = psp.zval
        c1 = psp.c[1]
        c2 = psp.c[2]
        c3 = psp.c[3]
        c4 = psp.c[4]
        v_of_0 = v_of_0 + 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15.0*c3 + 105.0*c4)  # to match QE
    end

    #v_of_0 = v_of_0 /Ham.pw.CellVolume #/prod(Ham.pw.Ns)
    v_of_0 = v_of_0 / Ham.pw.CellVolume * prod(Ham.pw.Ns) / Ham.pw.CellVolume

    println("v_of_0 = ", v_of_0)

    return v_of_0
end


function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
end

function init_Ham_CH4()
    # Atoms
    atoms = init_atoms_xyz(joinpath(DIR_STRUCTURES,"CH4.xyz"))
    atoms.LatVecs = gen_lattice_cubic(16.0)

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP,"C-q4.gth"),
                joinpath(DIR_PSP,"H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    return Ham
end

function init_Ham_GaAs_fcc( a::Float64, meshk::Array{Int64,1} )

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(a) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end


function main()
    
    #Ham = init_Ham_Si_fcc( 10.2631, [3,3,3] )
    #Ham = init_Ham_CH4()
    Ham = init_Ham_GaAs_fcc(10.6839444516, [3,3,3])

    v_of_0 = calc_v_of_0( Ham )


    # Si fcc
    #vref = -7.3656406697117621e-2 * prod(Ham.pw.Ns) / Ham.pw.CellVolume / 2
    #vref = -7.3656406697117621e-2 / 2

    # CH4
    #vref = -8.5397263096609246e-5 * prod(Ham.pw.Ns) / Ham.pw.CellVolume / 2

    # GaAs
    vref = 9.4502035797829315e-2 * prod(Ham.pw.Ns) / Ham.pw.CellVolume / 2

    println("vref = ", vref)
end

main()