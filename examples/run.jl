using LinearAlgebra
using Random
using Printf

using PWDFT

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

function driver()

    t1 = time_stamp("run.jl starts")
    
    Nargs = length(ARGS)
    @assert( Nargs >= 1 )

    str_prg = "include(\"" * ARGS[1] * "\")"
    
    println("str_prg = ", str_prg)
    eval( Meta.parse(str_prg) )

    if Nargs == 2
        str_prg = "@time main(method=\"" * ARGS[2] * "\")"
    else
        str_prg = "@time main()"
    end
    
    println("str_prg = ", str_prg)
    
    # Run the main function two times
    # Set random seed to a same state before calling the main function

    Random.seed!(1234)
    eval( Meta.parse(str_prg) )
    
    Random.seed!(1234)
    eval( Meta.parse(str_prg) )

    t2 = time_stamp("\nrun.jl ends")

    println()
    println("Elapsed time = ", float( (t2 - t1).value )/1000, " seconds")
end

driver()
