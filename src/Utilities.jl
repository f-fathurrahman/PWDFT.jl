include("all_pade_psp.jl")
include("all_pbe_psp.jl")

#
# simpson's rule integration. On input:
#   mesh = the number of grid points (should be odd)
#   func(i)= function to be integrated
#   rab(i) = r(i) * dr(i)/di * di
#
# For the logarithmic grid not including r=0 :
#   r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
#
# For the logarithmic grid including r=0 :
#   r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
#
# Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
# where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
#
# Adapted from QE
#
function integ_simpson(Npoints, f, rab)
    asum = 0.0
    r12 = 1.0/3.0
    f3 = f[1]*rab[1]*r12

    for i = 2:2:Npoints-1
        f1 = f3
        f2 = f[i]*rab[i]*r12
        f3 = f[i+1]*rab[i+1]*r12
        asum = asum + f1 + 4.0*f2 + f3
    end

    return asum
end



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
