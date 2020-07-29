import PWDFT: op_V_loc, op_V_Ps_loc

function op_V_loc( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )
    
    Nstates = size(psis.data[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    out = zeros_BlochWavefuncGamma(Ham)
    
    for ispin = 1:Nspin
        Ham.ispin = ispin
        out.data[ispin] = op_V_loc( Ham, psis.data[ispin] )
    end

    return out

end


function op_V_Ps_loc( Ham::Hamiltonian, psis::BlochWavefuncGamma )

    Nstates = size(psis.data[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    out = zeros_BlochWavefuncGamma(Ham)
    
    for ispin = 1:Nspin
        Ham.ispin = ispin
        out.data[ispin] = op_V_Ps_loc( Ham, psis.data[ispin] )
    end

    return out

end


# apply V_Ps_loc, Hartree, and XC potentials
function op_V_loc( Ham::HamiltonianGamma, psi::Array{ComplexF64,2} )
    ispin = Ham.ispin
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( Ham.pw, V_loc, psi )
end

# only apply V_Ps_loc
function op_V_Ps_loc( Ham::HamiltonianGamma, psi::Array{ComplexF64,2} )
    return op_V_loc( Ham.pw, Ham.potentials.Ps_loc, psi )
end

# apply general V_loc
# ik must be given to get information about
# mapping between psi in G-space to real space
function op_V_loc( pw::PWGridGamma, V_loc, psi::Array{ComplexF64,2} )

    Ns = pw.Ns
    CellVolume  = pw.CellVolume
    Npoints = prod(Ns)
    Nstates = size(psi,2)
    
    idx_gw2r = pw.gvecw.idx_gw2r
    idx_gw2rm = pw.gvecw.idx_gw2rm
    Ngw = pw.gvecw.Ngw

    ctmp = zeros(ComplexF64, Npoints) # reduce memory requirement
    Vpsi = zeros(ComplexF64, Ngw, Nstates)

    for ist in 1:Nstates
        ctmp .= 0.0 + im*0.0
        ctmp[1] = psi[1,ist]
        for igw in 2:Ngw
            idx = idx_gw2r[igw]
            ctmp[idx] = psi[igw,ist]
            idxm = idx_gw2rm[igw]
            ctmp[idxm] = conj(psi[igw,ist])
        end

        # get values of psi in real space grid
        G_to_R!(pw, ctmp)

        for ip = 1:Npoints
            ctmp[ip] = V_loc[ip]*ctmp[ip]
        end

        # Transform back to G-space
        R_to_G!(pw, ctmp)

        #Vpsi[:,ist] .= ctmp[idx_gw2r] # FIXME: use loops?
        Vpsi[1,ist] = ctmp[1]
        for igw in 2:Ngw
            idx = idx_gw2r[igw]
            Vpsi[igw,ist] = ctmp[idx]
        end
    end

    return Vpsi
end

#
# single-column version
#
function op_V_loc( Ham::HamiltonianGamma, psi::Array{ComplexF64,1} )
    ispin = Ham.ispin
    V_loc = @view Ham.potentials.Total[:,ispin]
    return op_V_loc( Ham.ik, Ham.pw, V_loc, psi )
end

function op_V_Ps_loc( Ham::HamiltonianGamma, psi::Array{ComplexF64,1} )
    return op_V_loc( Ham.pw, Ham.potentials.Ps_loc, psi )
end

function op_V_loc( pw::PWGridGamma, V_loc, psi::Array{ComplexF64,1} )
    Ns = pw.Ns
    CellVolume  = pw.CellVolume
    Npoints = prod(Ns)
    ctmp = zeros(ComplexF64, Npoints)
    
    idx_gw2r = pw.gvecw.idx_gw2r
    idx_gw2rm = pw.gvecw.idx_gw2rm
    Ngw = pw.gvecw.Ngw

    ctmp[1] = psi[1]
    for igw in 2:Ngw
        idx = idx_gw2r[igw]
        ctmp[idx] = psi[igw]
        idxm = idx_gw2rm[igw]
        ctmp[idxm] = conj(psi[igw])
    end

    # get values of psi in real space grid
    G_to_R!(pw, ctmp)

    for ip = 1:Npoints
        ctmp[ip,ist] = V_loc[ip]*ctmp[ip,ist]
    end

    R_to_G!(pw, ctmp)

    return ctmp[idx_gw2r]
end

