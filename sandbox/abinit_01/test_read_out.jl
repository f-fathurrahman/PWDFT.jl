using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
include(joinpath(DIR_PWDFT, "utilities", "ABINIT.jl"))

function main()
    
    Natoms = read_abinit_out_Natoms("TEMP_abinit/LOG_H2O")
    println("Natoms = ", Natoms)
    
    forces = read_abinit_forces("TEMP_abinit/LOG_H2O")
    println(forces)

    println("Pass here")
end

main()