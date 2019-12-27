# -----------------------------------------------------------------------------
# Generate "zeros" wavefunction
# -----------------------------------------------------------------------------

function zeros_CuBlochWavefunc( Ham::CuHamiltonian )
    return zeros_CuBlochWavefunc( Ham.pw, Ham.electrons )
end

function zeros_CuBlochWavefunc( pw::CuPWGrid, electrons::CuElectrons )
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return zeros_CuBlochWavefunc( pw, Nstates, Nspin )
end

function zeros_CuBlochWavefunc( pw::CuPWGrid, Nstates::Int64, Nspin::Int64 )
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    Nkspin = Nspin*Nkpt

    psiks = CuBlochWavefunc(undef,Nkspin)

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = CuArrays.zeros( ComplexF64, Ngw[ik], Nstates )
    end
    return psiks
end



# -----------------------------------------------------------------------------
# Generate random wavefunction
# -----------------------------------------------------------------------------

function rand_CuBlochWavefunc( Ham::CuHamiltonian )
    return rand_CuBlochWavefunc( Ham.pw, Ham.electrons.Nstates, Ham.electrons.Nspin )
end

function rand_CuBlochWavefunc( pw::CuPWGrid, electrons::CuElectrons )
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return rand_CuBlochWavefunc( pw, Nstates, Nspin )
end

function rand_CuBlochWavefunc( pw::CuPWGrid, Nstates::Int64, Nspin::Int64 )
    
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    Nkspin = Nspin*Nkpt

    psiks = CuBlochWavefunc(undef,Nkspin)

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = rand_CuWavefunc(Ngw[ik],Nstates)
    end

    return psiks
end

function rand_CuWavefunc( Nbasis, Nstates )
    return ortho_gram_schmidt( CuArrays.rand(ComplexF64,Nbasis,Nstates) )
end

#
# Adapted from pwscf
#=
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
        ortho_sqrt(psiks[ikspin])
        Hr = Hermitian(psiks[ikspin]' * op_H(Ham, psiks[ikspin]))
        evals, evecs = eigen(Hr)
        #Ham.electrons.ebands[:,ik] = evals
        psiks[ikspin] = psiks[ikspin]*evecs
    end
    end

    return psiks

end
=#
