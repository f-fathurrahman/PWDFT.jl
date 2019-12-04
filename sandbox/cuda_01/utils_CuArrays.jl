#
# XXX: Should only works for Hermitian matrix
#
import LinearAlgebra: eigen
function eigen( d_A::CuArray{ComplexF64,2,Nothing} )
    return CUSOLVER.heevd!('V', 'U', d_A)
end