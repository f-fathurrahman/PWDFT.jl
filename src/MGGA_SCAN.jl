using Libxc

function calc_KEdens!(
    ik::Int64,
    pw::PWGrid,
    psi::Array{ComplexF64,2},
    KEdens::Array{Float64,1}
)
    G = pw.gvec.G
    Ngw = pw.gvecw.Ngw[ik]
    idx_gw2r = pw.gvecw.idx_gw2r[ik]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    Npoints = prod(pw.Ns)
    Nstates = size(psi,2)

    ∇ψx = zeros(ComplexF64,Npoints) 
    ∇ψy = zeros(ComplexF64,Npoints)
    ∇ψz = zeros(ComplexF64,Npoints)
    
    for ist in 1:Nstates
        fill!(∇ψx, 0.0 + im*0.0)
        fill!(∇ψy, 0.0 + im*0.0)
        fill!(∇ψz, 0.0 + im*0.0)
        for igw = 1:Ngw
            ig = idx_gw2g[igw]
            ip = idx_gw2r[igw]
            ∇ψx[ip] = im*G[1,ig]*psi[igw,ist]
            ∇ψy[ip] = im*G[2,ig]*psi[igw,ist]
            ∇ψz[ip] = im*G[3,ig]*psi[igw,ist]
        end
        G_to_R!(pw, ∇ψx)
        G_to_R!(pw, ∇ψy)
        G_to_R!(pw, ∇ψz)
        #
        ∇ψx[:] = ∇ψx[:]*sqrt(Npoints/pw.CellVolume)*sqrt(Npoints)
        ∇ψy[:] = ∇ψy[:]*sqrt(Npoints/pw.CellVolume)*sqrt(Npoints)
        ∇ψz[:] = ∇ψz[:]*sqrt(Npoints/pw.CellVolume)*sqrt(Npoints)
        # FIXME: Need to add wk and Focc weight?
        for ip in 1:Npoints
            KEdens[ip] = KEdens[ip] + real( conj(∇ψx[ip])*∇ψx[ip] +
                conj(∇ψy[ip])*∇ψy[ip] + conj(∇ψz[ip])*∇ψz[ip] )
        end
    end
    #@views KEdens[:] = 0.5*KEdens[:]
    return 
end



function calc_epsxc_SCAN(
    xc_calc::LibxcXCCalculator,
    pw::PWGrid,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,1}
)

    @assert size(psiks,1) == 1
    @assert pw.gvecw.kpoints.Nkpt == 1

    FUNC_IDX = 263 # mgga x scan
    FUNC_IDC = 267 # mgga c scan

    #FUNC_IDX = 202 # mgga x tpss
    #FUNC_IDC = 231 # mgga c tpss

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # Need to symmetryize KEdens?
    KEdens = zeros(Float64,Npoints)
    calc_KEdens!(1, pw, psiks[1], KEdens)

    lapl = zeros(Npoints)

    # apply threshold
    #for ip in 1:Npoints
    #    #Rhoe[ip] = max(Rhoe[ip], 1e-12)
    #    Rhoe[ip] = abs(Rhoe[ip])
    #    gRhoe2[ip] = max(gRhoe2[ip], 1e-24)
    #    KEdens[ip] = max(KEdens[ip], 1e-12)
    #end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, FUNC_IDX, Nspin)
    Libxc_xc_mgga_exc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, FUNC_IDC, Nspin)
    Libxc_xc_mgga_exc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_Vxc_SCAN!(
    xc_calc::LibxcXCCalculator,
    pw::PWGrid,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,1},
    V_xc::Array{Float64,1}
)

    @assert size(psiks,1) == 1

    FUNC_IDX = 263 # mgga x scan
    FUNC_IDC = 267 # mgga c scan

    #FUNC_IDX = 202 # mgga x tpss
    #FUNC_IDC = 231 # mgga c tpss

    Npoints = size(Rhoe,1)
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # Need to symmetryize KEdens?
    KEdens = zeros(Npoints)
    calc_KEdens!(1, pw, psiks[1], KEdens)

    # Not used
    lapl = zeros(Npoints)
    Vlapl = zeros(Npoints)

    # apply threshold
    #for ip in 1:Npoints
    #    #Rhoe[ip] = max(Rhoe[ip], 1e-12)
    #    Rhoe[ip] = abs(Rhoe[ip])
    #    gRhoe2[ip] = max(gRhoe2[ip], 1e-24)
    #    KEdens[ip] = max(KEdens[ip], 1e-12)
    #end

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)

    Vtau_x = zeros(Float64,Npoints)
    Vtau_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, FUNC_IDX, Nspin)
    Libxc_xc_mgga_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, V_x, Vg_x, Vlapl, Vtau_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, FUNC_IDC, Nspin)
    Libxc_xc_mgga_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, V_c, Vg_c, Vlapl, Vtau_c)
    Libxc_xc_func_end(ptr)

    # gradient correction
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip = 1:Npoints
        hx[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[1,ip]
        hy[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[2,ip]
        hz[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip = 1:Npoints
        V_xc[ip] = V_x[ip] + V_c[ip] - 2.0*divh[ip]
        xc_calc.Vtau[ip,1] = Vtau_x[ip] + Vtau_c[ip]
    end

    return
end

