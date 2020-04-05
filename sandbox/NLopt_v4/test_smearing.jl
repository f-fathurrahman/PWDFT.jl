using PWDFT

include("smearing.jl")

function main()
    
    evals = [1.0, 2.0, 3.0, 4.0]
    E_f = 2.5
    kT = 0.1
    
    println( smear_fermi.( evals, E_f, kT ) )

    x = (E_f .- evals)/kT
    println( PWDFT.wgauss.(x) )
end

main()