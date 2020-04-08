using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_NLsolve.jl"))
include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

include("calc_energies_grad.jl")
include("create_Ham.jl")

function main()
    Random.seed!(1234)

    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    #Ham = create_Ham_Si_fcc()
    Ham = create_Ham_GaAs()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_ZnO()

    psiks = rand_BlochWavefunc( Ham )

    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)

    g = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )

    Ham.energies.NN = calc_E_NN( Ham.atoms )

    E0 = calc_energies_grad!( Ham, psiks, g, Kg )
    println("E0 = ", E0)

    println("dot_BlochWavefunc (g,g) = ", dot_BlochWavefunc(Ham.pw.gvecw.kpoints, g,g))
    println("dot(g,psiks) = ", 2*real(dot(g,psiks)))

    for i in 1:length(psiks)
        @printf("i = %d dot(psiks,g) = %18.10e\n", i, 2*real(dot(psiks[i],g[i])))
    end

    for i in 1:length(psiks)
        @printf("i = %d dot(g,g) = %18.10e\n", i, 2*real(dot(g[i],g[i])))
    end

    #d = -Kg
    d = -g
    #d = rand_BlochWavefunc( Ham )
    #for i in 1:length(psiks)
    #    ortho_sqrt!(d[i])
    #end

    println("Before constrain: ", 2.0*dot(d,d))
    constrain_search_dir!( d, psiks )
    println("After constrain: ", 2.0*dot(d,d))

    #for i in 1:length(psiks)
    #    c = real(dot(d[i],d[i]))
    #    d[i] = d[i]/sqrt(c)
    #end
    println("dot(psiks,psiks) = ", dot(psiks,psiks))
    println("dot(d,d) = ", dot(d,d))

    psic = zeros_BlochWavefunc( Ham )
    for α in 10.0 .^ range(1,stop=-10,step=-1)
        #
        dE = 2.0*real( dot(g, α*d) )
        #dE = dot_BlochWavefunc(Ham.pw.gvecw.kpoints, g, α*d) #*length(g)/prod(Ham.pw.gvecw.kpoints.mesh)
        #dE = real( dot(g, α*d) )
        @printf("α = %e, dE = %18.10e\n", α, dE)
        #
        psic = deepcopy( psiks )
        do_step!( psic, α, d )
        Etot = calc_energies_grad!( Ham, psic, g, Kg )
        ratio = (Etot - E0)/dE
        @printf("α = %e, ratio = %18.10e, diff = %e\n", α, ratio, abs(1-ratio))
    end

    println("Pass here")

end

main()