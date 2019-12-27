import PWDFT: op_V_loc, op_V_Ps_loc


function op_V_loc( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_CuBlochWavefunc(Ham)
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_loc( Ham, psiks[ikspin] )
    end
    return out
end

function op_V_Ps_loc( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_CuBlochWavefunc(Ham)
    
    for ispin = 1:Nspin
    for ik=1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_Ps_loc( Ham, psiks[ikspin] )
    end
    end
    return out
end

# apply V_Ps_loc, Hartree, and XC potentials
function op_V_loc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )
    ispin = Ham.ispin
    ik = Ham.ik
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( ik, Ham.pw, V_loc, psi )
end

# only apply V_Ps_loc
function op_V_Ps_loc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )
    return op_V_loc( Ham.ik, Ham.pw, Ham.potentials.Ps_loc, psi )
end



# apply general V_loc
# ik must be given to get information about
# mapping between psi in G-space to real space
function op_V_loc( ik::Int64, pw::CuPWGrid, V_loc, psi::CuArray{ComplexF64,2} )

    Ns = pw.Ns
    CellVolume  = pw.CellVolume
    Npoints = prod(Ns)
    Nstates = size(psi,2)

    ctmp = CuArrays.zeros(ComplexF64, Npoints, Nstates)
    idx = pw.gvecw.idx_gw2r[ik]
    Ngw_k = length(idx)

    Nthreads = 256
    Nblocks = ceil(Int64, Ngw_k/Nthreads)

    for ist in 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_copy_to_fft_grid_gw2r!( ist, idx, psi, ctmp )        
    end

    # get values of psi in real space grid
    G_to_R!( pw, ctmp )

    for ist in 1:Nstates
        ctmp[:,ist] = V_loc[:].*ctmp[:,ist]
    end

    R_to_G!( pw, ctmp )

    cVpsi = CuArrays.zeros(ComplexF64, Ngw_k, Nstates)

    for ist in 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_copy_from_fft_grid_gw2r!( ist, idx, ctmp, cVpsi )        
    end

    return cVpsi
end



#
# single-column version
#
function op_V_loc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,1} )
    ispin = Ham.ispin
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( Ham.ik, Ham.pw, V_loc, psi )
end


function op_V_Ps_loc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,1} )
    return op_V_loc( Ham.ik, Ham.pw, Ham.potentials.Ps_loc, psi )
end


function op_V_loc( ik::Int64, pw::CuPWGrid, V_loc, psi::CuArray{ComplexF64,1} )
    
    Ns = pw.Ns
    CellVolume = pw.CellVolume
    Npoints = prod(Ns)

    ctmp = CuArrays.zeros(ComplexF64, Npoints)
    idx = pw.gvecw.idx_gw2r[ik]
    Ngw_k = length(idx)

    Nthreads = 256
    Nblocks = ceil(Int64, Ngw_k/Nthreads)
    @cuda threads=Nthreads blocks=Nblocks kernel_copy_to_fft_grid_gw2r_1state!( idx, psi, ctmp )

    # get values of psi in real space grid
    G_to_R!( pw, ctmp )

    ctmp[:] = V_loc[:].*ctmp[:]

    R_to_G!( pw, ctmp )
    cVpsi = CuArrays.zeros(ComplexF64, Ngw_k)
    
    @cuda threads=Nthreads blocks=Nblocks kernel_copy_from_fft_grid_gw2r_1state!( idx, ctmp, cVpsi )    

    return cVpsi

end

