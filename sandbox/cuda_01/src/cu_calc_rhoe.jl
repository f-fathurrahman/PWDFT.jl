import PWDFT: calc_rhoe, calc_rhoe!

function calc_rhoe!( Ham::CuHamiltonian, psiks::CuBlochWavefunc, Rhoe::CuArray{Float64,2}; renormalize=true )

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nelectrons_true = Ham.electrons.Nelectrons

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    Npoints = prod(Ns)
    Nstates = size(psiks[1])[2]

    psiR = CuArrays.zeros(ComplexF64, Npoints, Nstates)

    # dont forget to zero out the Rhoe first
    Rhoe[:,:] .= 0.0

    Nthreads = 256

    for ispin = 1:Nspin, ik = 1:Nkpt

        ikspin = ik + (ispin - 1)*Nkpt

        psiR[:,:] .= 0.0 + im*0.0
        
        idx = pw.gvecw.idx_gw2r[ik]
        psi = psiks[ikspin]

        Nblocks = ceil( Int64, Ngw[ik]/Nthreads )

        for ist in 1:Nstates
            @cuda threads=Nthreads blocks=Nblocks kernel_copy_to_fft_grid_gw2r!( ist, idx, psi, psiR )
        end

        # Transform to real space
        G_to_R!( pw, psiR )

        # orthonormalization in real space
        #ortho_gram_schmidt!( psiR ) # is this needed or simply scaling ?
        #psiR[:] = sqrt(Npoints/CellVolume)*psiR[:]

        psiR[:,:] .= sqrt(Npoints/CellVolume)*sqrt(Npoints)*psiR[:,:] # by pass orthonormalization, only use scaling

        for ist in 1:Nstates
            w = wk[ik]*Focc[ist,ikspin]
            Rhoe[:,ispin] .= Rhoe[:,ispin] .+ w*real( conj(psiR[:,ist]) .* psiR[:,ist] )
        end
    end

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        Rhoe[:] = Nelectrons_true/integ_rho * Rhoe[:]
    end

    #
    # XXX This is rather difficult to parallelize
    #
    # Symmetrize Rhoe if needed
    #if Ham.sym_info.Nsyms > 1
    #    symmetrize_rhoe!( Ham.pw, Ham.sym_info, Ham.rhoe_symmetrizer, Rhoe )
    #end

    return
end

function calc_rhoe( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = CuArrays.zeros(Float64, Npoints, Nspin)
    calc_rhoe!( Ham, psiks, Rhoe )
    return Rhoe
end


