function op_V_loc( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_loc( Ham, psiks[ikspin] )
    end
    return out
end


# In-place, accumulated version
function op_V_loc!( Ham::Hamiltonian, psiks::BlochWavefunc, Hpsiks::BlochWavefunc )
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        op_V_loc!( Ham, psiks[i], Hpsiks[i] )
    end
    return
end


# In-place, accumated version
function op_V_loc!( Ham::Hamiltonian,
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
)

    pw = Ham.pw

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = pw.gvecw.kpoints.Nkpt
    i = ik + (ispin - 1)*Nkpt
    V_loc = Ham.potentials.Total
    
    CellVolume  = pw.CellVolume
    Npoints = prod(pw.Ns)
    Nstates = size(psi,2)
    idx = pw.gvecw.idx_gw2r[ik]
    Ngw_ik = pw.gvecw.Ngw[ik]
    #
    ctmp = zeros(ComplexF64, Npoints)
    #
    for ist in 1:Nstates
        #
        fill!(ctmp, 0.0 + im*0.0)
        for igw in 1:Ngw_ik
            ip = idx[igw]
            ctmp[ip] = psi[igw,ist]
        end
        #
        # get values of psi in real space grid
        #
        G_to_R!(pw, ctmp)
        #
        # Multiply in real space
        #
        for ip in 1:Npoints
            ctmp[ip] = V_loc[ip,ispin]*ctmp[ip]
        end
        #
        # Back to G-space
        #
        R_to_G!(pw, ctmp)
        #
        # Accumalate result in Hpsi
        #
        for igw in 1:Ngw_ik
            ip = idx[igw]
            Hpsi[igw,ist] = Hpsi[igw,ist] + ctmp[ip]
        end
    end
    return
end


function op_V_Ps_loc( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_Ps_loc( Ham, psiks[ikspin] )
    end
    return out
end

# apply V_Ps_loc, Hartree, and XC potentials
function op_V_loc( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,2} )
    ispin = Ham.ispin
    ik = Ham.ik
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( ik, Ham.pw, V_loc, psi )
end

# only apply V_Ps_loc
function op_V_Ps_loc( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,2} )
    return op_V_loc( Ham.ik, Ham.pw, Ham.potentials.Ps_loc, psi )
end

# apply general V_loc
# ik must be given to get information about
# mapping between psi in G-space to real space
function op_V_loc( ik::Int64, pw::PWGrid, V_loc,
    psi::AbstractArray{ComplexF64,2}
)

    Ns = pw.Ns
    CellVolume  = pw.CellVolume
    Npoints = prod(Ns)
    Nstates = size(psi,2)

    ctmp = zeros(ComplexF64, Npoints)
    Vpsi = zeros(ComplexF64, pw.gvecw.Ngw[ik], Nstates)
    idx = pw.gvecw.idx_gw2r[ik]
    for ist in 1:Nstates
        #ctmp .= 0.0 + im*0.0
        fill!(ctmp, 0.0 + im*0.0)
        #
        ctmp[idx] = psi[:,ist]
        # get values of psi in real space grid
        G_to_R!(pw, ctmp)
        # Multiply in real space
        for ip in 1:Npoints
            ctmp[ip] = V_loc[ip]*ctmp[ip]
        end
        # Back to G-space
        R_to_G!(pw, ctmp)
        #
        Vpsi[:,ist] = ctmp[idx]
    end

    return Vpsi
end

#
# single-column version
#
function op_V_loc( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,1} )
    ispin = Ham.ispin
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( Ham.ik, Ham.pw, V_loc, psi )
end

function op_V_Ps_loc( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,1} )
    return op_V_loc( Ham.ik, Ham.pw, Ham.potentials.Ps_loc, psi )
end

function op_V_loc( ik::Int64, pw::PWGrid, V_loc,
    psi::AbstractArray{ComplexF64,1}
)
    Ns = pw.Ns
    CellVolume  = pw.CellVolume
    Npoints = prod(Ns)

    ctmp = zeros(ComplexF64, Npoints)
    idx = pw.gvecw.idx_gw2r[ik]
    ctmp[idx] = psi[:]

    # get values of psi in real space grid
    ctmp = G_to_R( pw, ctmp )

    for ip in 1:Npoints
        ctmp[ip] = V_loc[ip]*ctmp[ip]
    end

    cVpsi = R_to_G( pw, ctmp )
    return cVpsi[idx]
end

