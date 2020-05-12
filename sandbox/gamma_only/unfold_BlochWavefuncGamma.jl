function unfold_BlochWavefuncGamma( pw::PWGridGamma, pw_full::PWGrid, psis::BlochWavefuncGamma )
    
    Nspin = length(psis)
    Ngw = size(psis.data[1],1)
    Nstates = size(psis.data[1],2)
    
    Npoints = prod(pw.Ns)
    ctmp = zeros(ComplexF64,Npoints)

    psiks = Vector{Array{ComplexF64,2}}(undef,Nspin)
    for ispin in 1:Nspin
        psiks[ispin] = zeros(ComplexF64, pw_full.gvecw.Ngw[ispin], Nstates)
    end

    for ispin in 1:Nspin
        
        psi = psis.data[ispin]
        
        idx_gw2r_full = pw_full.gvecw.idx_gw2r[ispin]

        for ist in 1:Nstates
            ctmp[1] = psi[1,ist]
            for igw in 2:Ngw
                ip = pw.gvecw.idx_gw2r[igw]
                ipm = pw.gvecw.idx_gw2rm[igw]
                ctmp[ip] = psi[igw,ist]
                ctmp[ipm] = conj(psi[igw,ist])
            end
            psiks[ispin][:,ist] = ctmp[idx_gw2r_full]
        end
    end

    return psiks

end