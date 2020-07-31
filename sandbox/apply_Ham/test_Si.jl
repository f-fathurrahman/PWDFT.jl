using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main( ; method="SCF" )

    Random.seed!(1234)

    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[8,8,8] )

    psiks = rand_BlochWavefunc(Ham)

    Hpsiks_ref = op_H(Ham, psiks)
    print("Using ref op_H: ")
    @time begin
        for i in 1:10
            Hpsiks_ref = op_H(Ham, psiks)
        end
    end

    Hpsiks = zeros_BlochWavefunc(Ham)
    op_H!(Ham, psiks, Hpsiks)
    print("Using in-place op_H: ")
    @time begin
        for i in 1:10
            op_H!(Ham, psiks, Hpsiks)
        end
    end

    Nkspin = length(psiks)
    for i in 1:Nkspin
        println("\ni = ", i)
        Ngw = length(Hpsiks[i])
        s1 = sum(Hpsiks[i])
        println("sum Hpsiks     = ", s1)
        s2 = sum(Hpsiks_ref[i])
        println("sum Hpsiks_ref = ", s2)
        println("avg diff = ", (s1 - s2)/Ngw, " (should be small)")
    end

    println()
    s1 = dot(psiks,Hpsiks)
    println("dot(psiks,Hpsiks)     = ", s1)
    s2 = dot(psiks,Hpsiks_ref)
    println("dot(psiks,Hpsiks_ref) = ", s2)
    println("diff = ", s1 - s2, " (should be small)")
end

main()