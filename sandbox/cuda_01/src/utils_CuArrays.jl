#
# XXX: Should only works for Hermitian matrix
#

function eigen_heevd_gpu( d_A::CuArray{ComplexF64,2} )
    d_A_c = copy(d_A)
    return CUSOLVER.heevd!('V', 'U', d_A_c)
end

function eigvals_heevd_gpu( d_A::CuArray{ComplexF64,2} )
    d_A_c = copy(d_A)
    return CUSOLVER.heevd!('N', 'U', d_A_c)
end