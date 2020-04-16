using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

function main()
    
    Natoms = read_pwscf_out_Natoms("TEMP_pwscf/LOG_H2O")
    println("Natoms = ", Natoms)
    forces = read_pwscf_forces("TEMP_pwscf/LOG_H2O")
    println(forces)

    println("Pass here")
end

main()