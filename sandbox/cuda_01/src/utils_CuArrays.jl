#
# XXX: Should only works for Hermitian matrix
#
import LinearAlgebra: eigen, eigvals

function eigen( d_A::CuArray{ComplexF64,2} )
    d_A_c = copy(d_A)
    return CUSOLVER.heevd!('V', 'U', d_A_c)
end

function eigvals( d_A::CuArray{ComplexF64,2} )
    d_A_c = copy(d_A)
    return CUSOLVER.heevd!('N', 'U', d_A_c)
end