using Printf
using LinearAlgebra

function test_vector()
    
    # the DC component (first element) is real
    
    v1 = [3.1, 3.0 + im*2, 4.1 + im*3, 3.0 - im*2, 4.1 - im*3]

    v1g = [3.1, 3.0 + im*2, 4.1 + im*3] # gamma only

    println("dot v1  = ", dot(v1,v1))
    println("dot v1g = ", dot(v1g,v1g))
    println("gamma only dot = ", 2*dot(v1g,v1g) - v1g[1]*v1g[1])
end

test_vector()