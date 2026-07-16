function _exx_qgrid_init!(
    pw, nqs, nq1, nq2, nq3,
    max_nk, temp_nkqs, temp_xkq, temp_index_ikq, dxk
)
    SMALL_Q = 1e-6 #XXX FIXME: Hardcoded

    Nkpt = pw.gvecw.kpoints.Nkpt

    new_ikq = zeros(Int64, max_nk)
    index_xkq = zeros(Int64, Nkpt, nqs)
    # nqs is computed with current nq1, nq2, nq3

    # should be returned
    nkqs = 0 

    sxk = zeros(Float64, 3)
    xk_cryst = zeros(Float64, 3)

    # define the q-mesh step-sizes
    dq1 = 1.0/nq1
    dq2 = 1.0/nq2
    dq3 = 1.0/nq3
    invRecVecs = inv(pw.RecVecs)
    for ik in 1:Nkpt
        # go to crystalline coordinates
        xk_cryst[:] = invRecVecs*pw.gvecw.kpoints.k[:,ik]
        #
        iq = 0
        for iq1 in 1:nq1, iq2 in 1:nq2, iq3 in 1:nq3
            sxk[1] = xk_cryst[1] + (iq1-1)*dq1
            sxk[2] = xk_cryst[2] + (iq2-1)*dq2
            sxk[3] = xk_cryst[3] + (iq3-1)*dq3
            iq += 1
            xk_not_found = true
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                        #@printf("%4d%4d%4d%18.10f%18.10f%18.10f\n", iq1, iq2, iq3, dxk...)
                        if new_ikq[ikq] == 0
                            nkqs += 1
                            temp_index_ikq[nkqs] = ikq
                            new_ikq[ikq] = nkqs
                        end
                        index_xkq[ik,iq] = new_ikq[ikq]
                    end # if
                end # if
            end # for ikq
            #
            if xk_not_found
                return nkqs, index_xkq
            end
        end
    end #
    return nkqs, index_xkq
end