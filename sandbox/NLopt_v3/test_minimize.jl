using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("calc_energies_grad.jl")
include("create_Ham.jl")
include("KS_solve_Emin_PCG_new.jl")

function main()
    Random.seed!(1234)
    
    Ham = create_Ham_H2()    
    psiks = rand_BlochWavefunc( Ham )
    
    KS_solve_Emin_PCG_new!( Ham, psiks )

end

main()
