function _exx_grid_check(
    pw, sym_info,
    nq1, nq2, nq3, index_xk, index_xkq, index_sym
)

    SMALL_Q = 1e-6 #XXX FIXME: Hardcoded

    s = sym_info.s
    Nkpt = pw.gvecw.kpoints.Nkpt
    LatVecs = pw.LatVecs

    sxk = zeros(Float64, 3)
    dxk = zeros(Float64, 3)
    xk_cryst = zeros(Float64, 3)
    xkk_cryst = zeros(Float64, 3)
  
    dq1 = 1.0/nq1
    dq2 = 1.0/nq2
    dq3 = 1.0/nq3

    #println("dq1=$dq1, dq2=$dq2, dq3=$dq3")
    invRecVecs = inv(pw.RecVecs)
    for ik in 1:Nkpt
        #println("\nBegin check ik = ", ik)
        xk_cryst[:] = invRecVecs*pw.gvecw.kpoints.k[:,ik] # to crystal coordinates
        #println("xk_cryst = ", xk_cryst)
        iq = 0
        for iq1 in 1:nq1, iq2 in 1:nq2, iq3 in 1:nq3
            sxk[1] = xk_cryst[1] + (iq1-1)*dq1
            sxk[2] = xk_cryst[2] + (iq2-1)*dq2
            sxk[3] = xk_cryst[3] + (iq3-1)*dq3
            iq += 1
            ikq = index_xkq[ik,iq]
            ikk = index_xk[ikq]
            isym = index_sym[ikq]
            #
            xkk_cryst[:] = LatVecs[1,:]*pw.gvecw.kpoints.k[1,ikk]/(2π) +
                           LatVecs[2,:]*pw.gvecw.kpoints.k[2,ikk]/(2π) +
                           LatVecs[3,:]*pw.gvecw.kpoints.k[3,ikk]/(2π)
            #println("xkk_cryst = ", xkk_cryst)
            if isym < 0
                xkk_cryst[:] = -xkk_cryst[:]
            end
            isym = abs(isym)
            dxk[:] = s[:,1,isym]*xkk_cryst[1] +
                     s[:,2,isym]*xkk_cryst[2] +
                     s[:,3,isym]*xkk_cryst[3] - sxk[:]
            #println("dxk before round = ", dxk)
            dxk[:] = dxk[:] - round.(Int64, dxk)
            #println("dxk after round = ", dxk)
            if !all( abs.(dxk) .<= SMALL_Q )
                println(ik,iq)
                println(ikq,ikk,isym)
                println(dxk)
                error("Something wrong in exx grid init")
            end
        end
    end
    return
end