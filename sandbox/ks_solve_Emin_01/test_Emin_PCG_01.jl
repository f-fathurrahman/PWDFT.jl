using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("KS_solve_Emin_PCG_01.jl")

function create_Ham_Si_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[4,4,4] )
end

function create_Ham_GaN()
    atoms = Atoms( xyz_string_frac=
        """
        4

        Ga   0.333333333333333   0.666666666666667   0.000000000000000 
        Ga   0.666666666666667   0.333333333333333   0.500000000000000 
         N   0.333333333333333   0.666666666666667   0.385000000000000 
         N   0.666666666666667   0.333333333333333   0.885000000000000
        """, in_bohr=true,
        LatVecs = gen_lattice_hexagonal( 3.18*ANG2BOHR, 5.166*ANG2BOHR ) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "N-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_GaAs()
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(10.6839444516) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 20.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[5,5,5] )
end

function create_Ham_ZnO()
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        4

        Zn      0.3333333   0.6666667   0.0000000
        Zn      0.6666667   0.3333333   0.5000000
        O       0.3333333   0.6666667   0.3450000
        O       0.6666667   0.3333333   0.8450000
        """, in_bohr=true,
        LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR ) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Zn-q2.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[4,4,4] )
end


function main( ; method="SCF" )

    #Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaN()
    #Ham = create_Ham_GaAs()
    Ham = create_Ham_ZnO()

    Random.seed!(1234)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "Emin_01"
        KS_solve_Emin_PCG_01!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=50 )

    else
        error( @sprintf("Unknown method %s", method) )
    end
    
end

@time main(method="Emin")
@time main(method="Emin_01")

#=
    Kinetic energy  =  3.21045499991613E+00
    Hartree energy  =  5.76194835366404E-01
    XC energy       = -2.41016569027563E+00
    Ewald energy    = -8.39792740071415E+00
    PspCore energy  = -2.94625629171302E-01
    Loc. psp. energy= -2.17613632002650E+00
    NL   psp  energy=  1.58119195731965E+00
    >>>>>>>>> Etotal= -7.91101324758539E+00

PWDFT.jl:
Kinetic    energy:       3.2107141913
Ps_loc     energy:      -2.1756827897
Ps_nloc    energy:       1.5804698683
Hartree    energy:       0.5829626011
XC         energy:      -2.4156830816
PspCore    energy:      -0.2946256268
-------------------------------------
Electronic energy:       0.4881551626
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.9097722381


QE:
one-elec   energy:       2.3207092600
Hartree    energy:       0.5764842900
XC         energy:      -2.4102748100
-------------------------------------
Electronic energy:       0.4869187400
NN         energy:      -8.3979274200
-------------------------------------
Total      energy:      -7.9110086800
=#
