function calc_rhoe!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,2};
    renormalize=true
)
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nelectrons_true = Ham.electrons.Nelectrons

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    idx_gw2r = pw.gvecw.idx_gw2r
    Npoints = prod(Ns)
    Nstates = size(psiks[1])[2]

    ctmp = zeros(ComplexF64, Npoints)

    # dont forget to zero out the Rhoe first
    fill!(Rhoe, 0.0)
    NptsPerSqrtVol = Npoints/sqrt(CellVolume)
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        #
        for ist in 1:Nstates
            #
            fill!(ctmp, 0.0 + im*0.0)
            #
            for igw in 1:Ngw[ik]
                ip = idx_gw2r[ik][igw]
                ctmp[ip] = psi[igw,ist]
            end
            # to real space
            G_to_R!(pw, ctmp)
            # Renormalize
            for ip in 1:Npoints
                ctmp[ip] *= NptsPerSqrtVol
            end
            w = wk[ik]*Focc[ist,ikspin]
            for ip in 1:Npoints
                # accumulate
                Rhoe[ip,ispin] += w*real( conj(ctmp[ip])*ctmp[ip] )
            end
        end
    end # ik, ispin

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        for i in 1:length(Rhoe)
            Rhoe[i] *= Nelectrons_true/integ_rho
        end
    end

    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham.pw, Ham.sym_info, Ham.rhoe_symmetrizer, Rhoe )
    end

    return
end

function calc_rhoe( Ham::Hamiltonian, psiks::BlochWavefunc )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe!( Ham, psiks, Rhoe )
    return Rhoe
end
