using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include(joinpath(DIR_PWDFT, "sandbox", "DATA_DeltaCodes", "read_pwscf_input.jl"))
const PWINPUT_DIR = joinpath(DIR_PWDFT, "sandbox", "DATA_DeltaCodes", "PWINPUT")

include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_NLsolve.jl"))

function init_Hamiltonian()
    atoms, meshk = read_pwscf_input(joinpath(PWINPUT_DIR, "Al.in"))
    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="VWN",
                        meshk=[8,8,8], extra_states=4 )
end


function main()
    Random.seed!(1234)
    Ham = init_Hamiltonian()
    KS_solve_SCF_NLsolve!( Ham, use_smearing=true, kT=0.01 )
end

main()
@time main()