function set_exx_buffer!(Ham, psiks)
    
    exx = Ham.exx
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    Ngw = Ham.pw.gvecw.Ngw

    idx_gw2r = exx.idx_gw2r
    planbw = exx.planbw
    Ns = exx.Ns
    nkqs = exx.nkqs
    index_xk = exx.index_xk
    index_sym = exx.index_sym
    rir = exx.rir

    Npoints_exx = prod(Ns)
    psic_exx = zeros(ComplexF64, Npoints_exx)
    temppsic = zeros(ComplexF64, Npoints_exx)

    for ik in 1:Nkpt
        Ngwk = Ngw[ik]
        psi = psiks[ik]
        for ist in 1:Nstates # XXX should loop over active states for EXX
            fill!(temppsic, 0.0 + im*0.0)
            for igw in 1:Ngwk
                ip = idx_gw2r[ik][igw]
                temppsic[ip] = psi[igw,ist]
            end
            do_fft!(planbw, Ns, temppsic)
            #
            for ikq in 1:nkqs
                if index_xk[ikq] != ik
                    continue
                end
                #
                isym = abs(index_sym[ikq])
                #
                for ir in 1:Npoints_exx
                    psic_exx[ir] = temppsic[rir[ir,isym]]
                end
                #
                for ir in 1:Npoints_exx
                    if index_sym[ikq] < 0
                        psic_exx[ir] = conj(psic_exx[ir])
                    end
                    exx.buff[ir,ist,ikq] = psic_exx[ir]
                end
            end # ikq
        end # ist
    end # ik
    return
end