import InteractiveUtils
using Printf
using Random
using BenchmarkTools
using PWDFT

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_TIMING = joinpath(DIR_PWDFT, "timing")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

import Dates
function time_stamp(message::String)
    t1 = Dates.now()
    print(message, " : ")
    print(Dates.dayname(t1), ", ")
    print(Dates.day(t1), " ")
    print(Dates.monthname(t1), " ")
    print(Dates.year(t1), " ")
    print(Dates.hour(t1), ":")
    print(Dates.minute(t1), ":")
    print(Dates.second(t1))
    print("\n")
    return t1
end

println()
time_stamp("timing start:")
println()
println("Running on ", gethostname())
println()
InteractiveUtils.versioninfo()

const Nargs = length(ARGS)

if Nargs == 0
    const FILELIST = [
        "time_gen_lattice.jl",
        "time_Atoms.jl",
        "time_PWGrid.jl",
        "time_SymmetryInfo.jl",
        "time_KPoints.jl",
        "time_calc_rhoe.jl",
        "time_Hamiltonian_H.jl",
        "time_Hamiltonian_Pt.jl",
        "time_op_kpt.jl",
        "time_op_kpt_PBE.jl",
    ]    
    for fil in FILELIST
        include(joinpath(DIR_TIMING, fil))
    end
else
    include(ARGS[1])
end
