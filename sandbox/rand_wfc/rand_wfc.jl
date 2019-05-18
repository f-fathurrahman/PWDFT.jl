function rand_wfc( Ham::Hamiltonian )

    Nstates = Ham.electrons.Nstates
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw

    psiks = zeros_BlochWavefunc(Ham)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        for ist = 1:Nstates
            for igk = 1:Ngw[ik]
                ig = idx_gw2g[ik][igk]
                rr = rand()
                arg = 2*pi*rand()
                num = rr*cos(arg) + im*rr*sin(arg)
                Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                denum = Gw2 + 1.0
                psiks[ikspin][igk,ist] = num/denum
            end
        end

        Ham.ispin = ispin
        Ham.ik = ik

        Hr = Hermitian(psiks[ikspin]' * op_H(Ham, psiks[ikspin]))
        evals, evecs = eigen(Hr)
        #Ham.electrons.ebands[:,ik] = evals
        psiks[ikspin] = evecs[:,:]

        ortho_check(psiks[ikspin])
    end
    end

    return psiks

end
 