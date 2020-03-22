using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("calc_energies_grad.jl")
include("create_Ham.jl")

function main_test()
    Random.seed!(1234)
    Ham = create_Ham_H2()
    psiks = rand_BlochWavefunc( Ham )
    
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    
    energies = calc_energies( Ham, psiks )
    println("sum energies = ", sum(energies))
    
    Ham.ik = 1
    Ham.ispin = 1
    g_1 = calc_grad( Ham, psiks[1] )
    println("sum g_1 = ", sum(g_1))
    println("sum Kprec(g_1) = ", sum(Kprec(1, Ham.pw, g_1)))
end

function main()
    Random.seed!(1234)
    
    Ham = create_Ham_H2()
    
    psiks = rand_BlochWavefunc( Ham )
    
    g = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )

    Etot = calc_energies_grad!( Ham, psiks, g, Kg )
    println("Etot = ", Etot)

    println("sum g[1] = ", sum(g[1]))
    println("sum Kg[1] = ", sum(Kg[1]))
end

main()
main_test()