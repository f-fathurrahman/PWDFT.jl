include("all_pade_psp.jl")
include("all_pbe_psp.jl")

function inv_m3x3( A::Array{Float64,2} )

    @assert( (size(A)[1] == 3) && (size(A)[2] == 3) )

    # inverse determinant of the matrix
    detinv = 1.0/( A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2]
                 - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1]
                 + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1] )

    B = zeros(Float64,3,3)
    
    # inverse
    B[1,1] =  detinv * (A[2,2]*A[3,3] - A[2,3]*A[3,2])
    B[2,1] = -detinv * (A[2,1]*A[3,3] - A[2,3]*A[3,1])
    B[3,1] =  detinv * (A[2,1]*A[3,2] - A[2,2]*A[3,1])

    B[1,2] = -detinv * (A[1,2]*A[3,3] - A[1,3]*A[3,2])
    B[2,2] =  detinv * (A[1,1]*A[3,3] - A[1,3]*A[3,1])
    B[3,2] = -detinv * (A[1,1]*A[3,2] - A[1,2]*A[3,1])

    B[1,3] =  detinv * (A[1,2]*A[2,3] - A[1,3]*A[2,2])
    B[2,3] = -detinv * (A[1,1]*A[2,3] - A[1,3]*A[2,1])
    B[3,3] =  detinv * (A[1,1]*A[2,2] - A[1,2]*A[2,1])

    return B
end


function invTrans_m3x3( A::Array{Float64,2} )

    @assert( (size(A)[1] == 3) && (size(A)[2] == 3) )

    # inverse determinant of the matrix
    detinv = 1.0/( A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2]
                 - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1]
                 + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1] )

    B = zeros(Float64,3,3)
    
    # inverse and transpose
    B[1,1] =  detinv * (A[2,2]*A[3,3] - A[2,3]*A[3,2])
    B[1,2] = -detinv * (A[2,1]*A[3,3] - A[2,3]*A[3,1])
    B[1,3] =  detinv * (A[2,1]*A[3,2] - A[2,2]*A[3,1])

    B[2,1] = -detinv * (A[1,2]*A[3,3] - A[1,3]*A[3,2])
    B[2,2] =  detinv * (A[1,1]*A[3,3] - A[1,3]*A[3,1])
    B[2,3] = -detinv * (A[1,1]*A[3,2] - A[1,2]*A[3,1])

    B[3,1] =  detinv * (A[1,2]*A[2,3] - A[1,3]*A[2,2])
    B[3,2] = -detinv * (A[1,1]*A[2,3] - A[1,3]*A[2,1])
    B[3,3] =  detinv * (A[1,1]*A[2,2] - A[1,2]*A[2,1])

    return B
end


function transpose_m3x3( A::Array{Float64,2} )
    
    @assert( (size(A)[1] == 3) && (size(A)[2] == 3) )

    B = zeros(Float64,3,3)

    B[1,1] = A[1,1]
    B[1,2] = A[2,1]
    B[1,3] = A[3,1]

    B[2,1] = A[1,2]
    B[2,2] = A[2,2]
    B[2,3] = A[3,2]

    B[3,1] = A[1,3]
    B[3,2] = A[2,3]
    B[3,3] = A[3,3]
    
    return B
end


function print_matrix( A::Array{Float64,2} )
    Nrows = size(A)[1]
    Ncols = size(A)[2]
    for ir = 1:Nrows
        for ic = 1:Ncols
            @printf("%6.3f ", A[ir,ic])
        end
        @printf("\n")
    end
end

function print_matrix( A::Array{ComplexF64,2} )
    Nrows = size(A)[1]
    Ncols = size(A)[2]
    for ir = 1:Nrows
        for ic = 1:Ncols
            @printf("(%6.3f + %6.3fim) ", real(A[ir,ic]), imag(A[ir,ic]))
        end
        @printf("\n")
    end
end
