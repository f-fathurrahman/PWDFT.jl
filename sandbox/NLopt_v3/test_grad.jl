using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("calc_energies_grad.jl")

function create_Ham_GaAs_v1()
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(10.6839444516) )
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_GaAs_v2()
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(10.6839444516) )
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], use_symmetry=false, time_reversal=false )
end


function main()
    Ham1 = create_Ham_GaAs_v1()
    Ham2 = create_Ham_GaAs_v2()

    println(Ham1.pw.gvecw.kpoints)

    println(Ham2.pw.gvecw.kpoints)

    Nstates = Ham1.electrons.Nstates
    Nelectrons = Ham1.electrons.Nelectrons

    psiks1 = zeros_BlochWavefunc(Ham1)
    psiks2 = zeros_BlochWavefunc(Ham2)

    for i in 1:length(psiks1)
        for ist in 1:Nstates
            psiks1[i][:,ist] .= 10.0 + ist
        end
        ortho_sqrt!(psiks1[i])
        println("size(psiks1[i] = ", size(psiks1[i]))
        @printf("%2d %18.10f\n", i, real(sum(psiks1[i])))
    end

    g1 = zeros_BlochWavefunc(Ham1)
    Kg1 = zeros_BlochWavefunc(Ham1)
    Rhoe1 = calc_rhoe(Ham1, psiks1)
    update!( Ham1, Rhoe1 )
    Etot1 = calc_energies_grad!( Ham1, psiks1, g1, Kg1 )

    println("Etot1 = ", Etot1)
    println("dot(g1,g1) = ", dot(g1,g1))
    #wk_full = [1.0, 4.0, 1.0, 4.0, 2.0, 2.0, 2.0, 2.0]
    #ss = 0.0 + im*0.0
    for i in 1:length(psiks1)
        xx = dot(g1[i],g1[i]) #/wk_full[i]
        #ss = ss + xx
        println("i = ", i, " dot(g1,g1) = ", xx)
    end
    #println("2*ss = ", 2*ss)


    #i = 1
    #println("i = ", i, " dot(g1,g1)   = ", dot(g1[i],g1[i])/1)
    #i = 2
    #println("i = ", i, " dot(g1,g1)/4 = ", dot(g1[i],g1[i])/4)
    #i = 3
    #println("i = ", i, " dot(g1,g1)/3 = ", dot(g1[i],g1[i])/3)
    ##end
    #println("sum dot(g1,g1)   = ", dot(g1[1],g1[1]) + dot(g1[2],g1[2])/4 + dot(g1[3],g1[3])/3)


    # Second Hamiltonian

    for i in 1:length(psiks2)
        for ist in 1:Nstates
            psiks2[i][:,ist] .= 10.0 + ist
        end
        ortho_sqrt!(psiks2[i])
        #@printf("%2d %18.10f\n", i, real(sum(psiks2[i])))
    end

    g2 = zeros_BlochWavefunc(Ham2)
    Kg2 = zeros_BlochWavefunc(Ham2)
    Rhoe2 = calc_rhoe(Ham2, psiks2)
    update!( Ham2, Rhoe2 )
    Etot2 = calc_energies_grad!( Ham2, psiks2, g2, Kg2 )

    println("Etot2 = ", Etot2)
    println("dot(g2,g2) = ", dot(g2,g2))
    ss = 0.0 + im*0.0
    for i in 1:length(psiks2)
        xx = dot(g2[i],g2[i])
        ss = ss + xx
        println("i = ", i, " dot(g2,g2) = ", dot(g2[i],g2[i]))
    end
    println("2*ss = ", 2*ss)

    println("dot_BlochWavefunc(g1,g1) = ", dot_BlochWavefunc(Ham1.pw.gvecw.kpoints, g1,g1))
    println("dot_BlochWavefunc(g2,g2) = ", dot_BlochWavefunc(Ham2.pw.gvecw.kpoints, g2,g2))
    println("Pass here")
end

main()

