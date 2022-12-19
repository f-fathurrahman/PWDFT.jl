# This function can be merged with dense_to_smooth,
# however we make the operation explicit
function smooth_to_dense!(
    pw::PWGrid, v_in, v_out
)

    # Immediate return if not using dual grid
    if !pw.using_dual_grid
        v_out[:] = v_in[:]
        return
    end

    NptsDense = prod(pw.Ns)
    NptsSmooth = prod(pw.Nss) # not used?

    aux_in = zeros(ComplexF64, NptsSmooth)

    # Copy all data to aux_in
    aux_in[1:NptsSmooth] .= v_in[1:NptsSmooth]

    R_to_G!(pw, aux_in, smooth=true) # using smooth grid

    aux_out = zeros(ComplexF64, NptsDense)
    Ngs = pw.gvecs.Ng
    # Copy available plane wave coefs
    # For dense grid additional coefs are set to zero
    for ig in 1:Ngs
        ip_out = pw.gvec.idx_g2r[ig] # dense
        ip_in  = pw.gvecs.idx_g2r[ig] # smooth
        aux_out[ip_out] = aux_in[ip_in]
    end
    
    G_to_R!(pw, aux_out) # using dense grid

    # FIXME: Manual normalization
    aux_out *= (NptsDense/NptsSmooth)  # XXX This is important!

    v_out[1:NptsDense] = real(aux_out[1:NptsDense])

    return
end