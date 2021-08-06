function op_Vtau( Ham::Hamiltonian, psiks::BlochWavefunc )
    Vpsiks = zeros_BlochWavefunc(Ham)
    if Ham.xcfunc != "SCAN"
        return Vpsiks
    end
    Nkpt = size(psiks,1)
    # FIXME: spinpol is not yet implemented
    for ik in 1:Nkpt
        Ham.ik = ik
        op_Vtau!(Ham, psiks[ik], Vpsiks[ik])
    end
    return Vpsiks
end

function op_Vtau!( Ham::Hamiltonian, psiks::BlochWavefunc, Vpsiks::BlochWavefunc )
    if Ham.xcfunc != "SCAN"
        return
    end
    # FIXME: spinpol is not yet implemented
    Nkpt = size(psiks,1)
    for ik in 1:Nkpt
        Ham.ik = ik
        op_Vtau!(Ham, psiks[ik], Vpsiks[ik])
    end
    return Vpsiks
end

function op_Vtau( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,2} )
    Vpsi = zeros(ComplexF64,size(psi))
    op_Vtau!(Ham, psi, Vpsi)
    return Vpsi
end

# operator Vtau
function op_Vtau!( Ham::Hamiltonian,
    psi::AbstractArray{ComplexF64,2},
    Vpsi::AbstractArray{ComplexF64,2}
)
    if Ham.xcfunc != "SCAN"
        return
    end

    ik = Ham.ik

    pw = Ham.pw
    G = pw.gvec.G
    Ngw = pw.gvecw.Ngw[ik]
    idx_gw2r = pw.gvecw.idx_gw2r[ik]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    Npoints = prod(pw.Ns)
    Nstates = size(psi,2)
    Vtau = Ham.xc_calc.Vtau

    # use linear indexing
    Cx = zeros(ComplexF64, pw.Ns)
    Cy = zeros(ComplexF64, pw.Ns)
    Cz = zeros(ComplexF64, pw.Ns)
    
    for ist in 1:Nstates
        fill!(Cx, 0.0 + im*0.0)
        fill!(Cy, 0.0 + im*0.0)
        fill!(Cz, 0.0 + im*0.0)
        for igw = 1:Ngw
            ig = idx_gw2g[igw]
            ip = idx_gw2r[igw]
            Cx[ip] = im*G[1,ig]*psi[igw,ist]
            Cy[ip] = im*G[2,ig]*psi[igw,ist]
            Cz[ip] = im*G[3,ig]*psi[igw,ist]
        end
        G_to_R!(pw, Cx)
        G_to_R!(pw, Cy)
        G_to_R!(pw, Cz)
        # Multiply with Vtau in real space
        for ip in 1:Npoints
            Cx[ip] = Vtau[ip]*Cx[ip]
            Cy[ip] = Vtau[ip]*Cy[ip]
            Cz[ip] = Vtau[ip]*Cz[ip]
        end
        # to G-space
        R_to_G!(pw, Cx)
        R_to_G!(pw, Cy)
        R_to_G!(pw, Cz)
        # calculate the ∇⋅ in G-space
        # accumulate in Vpsi (with -0.5 factor)
        for igw in 1:Ngw
            ig = idx_gw2g[igw]
            ip = idx_gw2r[igw]
            Vpsi[igw,ist] = Vpsi[igw,ist] - 0.5*im*( G[1,ig]*Cx[ip] + G[2,ig]*Cy[ip] + G[3,ig]*Cz[ip] )
        end
    end
    return 
end