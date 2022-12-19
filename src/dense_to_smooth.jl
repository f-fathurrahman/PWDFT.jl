function dense_to_smooth!(
    pw::PWGrid, v_in, v_out
)

    # Immediate return if not using dual grid
    if !pw.using_dual_grid
        v_out[:] = v_in[:]
        return
    end

    NptsDense = prod(pw.Ns)
    NptsSmooth = prod(pw.Nss) # not used?

    aux_in = zeros(ComplexF64, NptsDense)

    # Copy all data to aux_in
    aux_in[1:NptsDense] .= v_in[1:NptsDense]

    R_to_G!(pw, aux_in) # using dense grid

    aux_out = zeros(ComplexF64, NptsSmooth)
    Ngs = pw.gvecs.Ng
    for ig in 1:Ngs
        ip_out = pw.gvecs.idx_g2r[ig] # smooth
        ip_in  = pw.gvec.idx_g2r[ig] # dense
        aux_out[ip_out] = aux_in[ip_in]
    end
    
    G_to_R!(pw, aux_out, smooth=true)
    # FIXME: Manual normalization
    aux_out *= (NptsSmooth/NptsDense)  # XXX This is important!

    v_out[1:NptsSmooth] = real(aux_out[1:NptsSmooth])

    return
end