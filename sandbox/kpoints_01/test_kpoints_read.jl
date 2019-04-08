using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()

    LATCONST = 10.2631

    atoms = Atoms(xyz_string_frac=
        """
        1

        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(LATCONST))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]

    kpoints = PWDFT.kpoints_from_string(atoms, """
    4
    0.000000000000   0.00000000000   0.00000000000   0.074074074074
    0.000000000000   0.00000000000   0.33333333333   0.592592592592
    0.000000000000   0.33333333333   0.33333333333   0.444444444444
    0.000000000000   0.33333333333  -0.33333333333   0.888888888889
    """)

    println(kpoints)


    #parse_kpts_string("""
    #4
    #0.0000000   0.0000000   0.0000000   0.0740741
    #0.0000000   0.0000000   0.3333333   0.5925926
    #0.0000000   0.3333333   0.3333333   0.4444444
    #0.0000000   0.3333333  -0.3333333   0.8888889
    #""")

    #Nkpt = 4
    #kpts = zeros(3,Nkpt)
    #wk = zeros(Nkpt)

    # k(    1) = (   0.0000000   0.0000000   0.0000000), wk =   0.0740741
    # k(    2) = (   0.0000000   0.0000000   0.3333333), wk =   0.5925926
    # k(    3) = (   0.0000000   0.3333333   0.3333333), wk =   0.4444444
    # k(    4) = (   0.0000000   0.3333333  -0.3333333), wk =   0.8888889

    #RecVecs = 2*pi*inv(Matrix(atoms.LatVecs'))
    ## kpts in crystal coord
    #kpts[:,1] = [0.0000000, 0.0000000,  0.0000000]   
    #kpts[:,2] = [0.0000000, 0.0000000,  0.3333333]   
    #kpts[:,3] = [0.0000000, 0.3333333,  0.3333333]   
    #kpts[:,4] = [0.0000000, 0.3333333, -0.3333333]   

    #wk = [0.0740741, 0.5925926, 0.4444444, 0.8888889]
    #ss = sum(wk)
    #wk = wk/ss  # normalize to 1

    #kpoints = KPoints(Nkpt, (0,0,0), RecVecs*kpts, wk, RecVecs)

    #println(kpoints)
end


#=
     number of k points=     4
                       cart. coord. in units 2pi/alat
        k(    1) = (   0.0000000   0.0000000   0.0000000), wk =   0.0740741
        k(    2) = (  -0.2357023   0.2357023  -0.2357023), wk =   0.5925926
        k(    3) = (   0.0000000   0.4714045   0.0000000), wk =   0.4444444
        k(    4) = (   0.4714045  -0.0000000   0.4714045), wk =   0.8888889

                       cryst. coord.
        k(    1) = (   0.0000000   0.0000000   0.0000000), wk =   0.0740741
        k(    2) = (   0.0000000   0.0000000   0.3333333), wk =   0.5925926
        k(    3) = (   0.0000000   0.3333333   0.3333333), wk =   0.4444444
        k(    4) = (   0.0000000   0.3333333  -0.3333333), wk =   0.8888889
=#


main()

