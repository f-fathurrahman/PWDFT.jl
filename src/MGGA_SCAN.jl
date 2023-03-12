using Libxc

function calc_KEdens!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    KEdens::Array{Float64,2}
)

    pw = Ham.pw
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    Focc = Ham.electrons.Focc
    k = pw.gvecw.kpoints.k

    G = pw.gvec.G
    Ngw = pw.gvecw.Ngw
    idx_gw2r = pw.gvecw.idx_gw2r
    idx_gw2g = pw.gvecw.idx_gw2g
    Npoints = prod(pw.Ns)
    Nstates = size(psiks[1],2)

    Nspin = size(KEdens, 2)
    @assert Nspin == 1

    fill!(KEdens, 0.0)

    ∇ψx = zeros(ComplexF64, Npoints) 
    ∇ψy = zeros(ComplexF64, Npoints)
    ∇ψz = zeros(ComplexF64, Npoints)
    
    # XXX: Assumption: Nspin == 1
    ispin = 1
    for ik in 1:Nkpt

        ikspin = ik + (ispin-1)*Nkpt
        psi = psiks[ikspin]

        for ist in 1:Nstates
            #
            fill!(∇ψx, 0.0 + im*0.0)
            fill!(∇ψy, 0.0 + im*0.0)
            fill!(∇ψz, 0.0 + im*0.0)
            #
            for igw in 1:Ngw[ik]
                ig = idx_gw2g[ik][igw]
                ip = idx_gw2r[ik][igw]
                ∇ψx[ip] = im*(G[1,ig] + k[1,ik])*psi[igw,ist]
                ∇ψy[ip] = im*(G[2,ig] + k[2,ik])*psi[igw,ist]
                ∇ψz[ip] = im*(G[3,ig] + k[3,ik])*psi[igw,ist]
            end
            #
            G_to_R!(pw, ∇ψx)
            G_to_R!(pw, ∇ψy)
            G_to_R!(pw, ∇ψz)
            #
            # Rescale
            scal = Npoints/sqrt(pw.CellVolume)
            @views ∇ψx[:] = ∇ψx[:]*scal
            @views ∇ψy[:] = ∇ψy[:]*scal
            @views ∇ψz[:] = ∇ψz[:]*scal
            # FIXME: Need to add wk and Focc weight?
            for ip in 1:Npoints
                KEdens[ip,ispin] += 0.5*Focc[ist]*wk[ik]*real(
                    conj(∇ψx[ip])*∇ψx[ip] + conj(∇ψy[ip])*∇ψy[ip] + conj(∇ψz[ip])*∇ψz[ip]
                )
            end
        end
    end

    # Symmetrize Rhoe if needed
    #println("Before symmetrize_rhoe: sum(KEdens) = ", sum(KEdens))
    #println("KEdens[1,1] = ", KEdens[1,1])
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( pw, Ham.sym_info, Ham.rhoe_symmetrizer, KEdens )
    end
    #println("After symmetrize_rhoe: sum(KEdens) = ", sum(KEdens))
    #println("KEdens[1,1] = ", KEdens[1,1])

    return 
end


# Evaluate both epsxc and Vxc for a metaGGA function (for now, only SCAN)
#
# Besides epsxc and V_xc, xc_calc.Vtau is also modified
# V_xc contains local and gradient correction
# xc_calc.Vtau contains the kinetic energy density correction
function calc_epsxc_Vxc_SCAN!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,1},
    epsxc::Array{Float64,1},
    V_xc::Array{Float64,1}
)

    fill!(epsxc, 0.0)
    fill!(V_xc, 0.0)

    xc_calc = Ham.xc_calc
    pw = Ham.pw

    # These are hardcoded for now
    # They should be taken from xc_calc.x_id and xc_calc.c_id
    FUNC_ID_X = 263 # mgga x scan
    FUNC_ID_C = 267 # mgga c scan

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2: |∇ρ|^2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip in 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]^2 + gRhoe[2,ip]^2 + gRhoe[3,ip]^2
    end

    # Calculate kinetic energy density here
    KEdens_ = zeros(Float64, Npoints, Nspin)
    calc_KEdens!(Ham, psiks, KEdens_)
    KEdens = reshape(KEdens_, Npoints) # need to be reshaped, for now

    lapl = zeros(Float64, Npoints) # Laplacian, not used

    # No threshold is applied


    # Memory for outputs
    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)

    Vtau_x = zeros(Float64,Npoints)
    Vtau_c = zeros(Float64,Npoints)

    Vlapl = zeros(Float64, Npoints) # Laplacian, not used

    #
    # Libxc is called here
    #
    ptr = Libxc_xc_func_alloc()
    #
    # exchange part
    #
    Libxc_xc_func_init(ptr, FUNC_ID_X, Nspin)
    @views Libxc_xc_mgga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_x, V_x, Vg_x, Vlapl, Vtau_x)
    Libxc_xc_func_end(ptr)
    #
    # correlation part
    #
    Libxc_xc_func_init(ptr, FUNC_ID_C, Nspin)
    @views Libxc_xc_mgga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_c, V_c, Vg_c, Vlapl, Vtau_c)
    Libxc_xc_func_end(ptr)
    #
    Libxc_xc_func_free(ptr)

    #
    # Energy
    #
    @views epsxc[:] .= eps_x[:] + eps_c[:]

    #
    # Potentials
    #

    # gradient correction
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip in 1:Npoints
        hx[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[1,ip]
        hy[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[2,ip]
        hz[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip in 1:Npoints
        V_xc[ip] = V_x[ip] + V_c[ip] - 2.0*divh[ip]
        xc_calc.Vtau[ip,1] = Vtau_x[ip] + Vtau_c[ip]
    end

    return

end





function calc_epsxc_SCAN(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,1}
)

    xc_calc = Ham.xc_calc
    pw = Ham.pw

    FUNC_ID_X = 263 # mgga x scan
    FUNC_ID_C = 267 # mgga c scan

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip in 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # Calculate kinetic energy density here
    KEdens_ = zeros(Float64, Npoints, Nspin)
    calc_KEdens!(Ham, psiks, KEdens_)
    KEdens = reshape(KEdens_, Npoints)

    lapl = zeros(Npoints)

#=
    # apply threshold
    for ip in 1:Npoints
        Rhoe[ip] = max(Rhoe[ip], 1e-12)
        #Rhoe[ip] = abs(Rhoe[ip])
        gRhoe2[ip] = max(gRhoe2[ip], 1e-24)
        KEdens[ip] = max(KEdens[ip], 1e-12)
    end
=#

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, FUNC_ID_X, Nspin)
    @views Libxc_xc_mgga_exc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, FUNC_ID_C, Nspin)
    @views Libxc_xc_mgga_exc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_Vxc_SCAN!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,1},
    V_xc::Array{Float64,1}
)

    xc_calc = Ham.xc_calc
    pw = Ham.pw

    FUNC_ID_X = 263 # mgga x scan
    FUNC_ID_C = 267 # mgga c scan

    Npoints = size(Rhoe,1)
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip in 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # Calculate kinetic energy density here
    KEdens_ = zeros(Npoints, Nspin)
    calc_KEdens!(Ham, psiks, KEdens_)
    KEdens = reshape(KEdens_, Npoints)

    # Not used
    lapl = zeros(Npoints)
    Vlapl = zeros(Npoints)

#=
    # apply threshold
    for ip in 1:Npoints
        Rhoe[ip] = max(Rhoe[ip], 1e-12)
    #    Rhoe[ip] = abs(Rhoe[ip])
        gRhoe2[ip] = max(gRhoe2[ip], 1e-24)
        KEdens[ip] = max(KEdens[ip], 1e-12)
    end
=#

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)

    Vtau_x = zeros(Float64,Npoints)
    Vtau_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, FUNC_ID_X, Nspin)
    Libxc_xc_mgga_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, V_x, Vg_x, Vlapl, Vtau_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, FUNC_ID_C, Nspin)
    Libxc_xc_mgga_vxc!(ptr, Npoints, Rhoe, gRhoe2, lapl, KEdens, V_c, Vg_c, Vlapl, Vtau_c)
    Libxc_xc_func_end(ptr)

    # gradient correction
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip in 1:Npoints
        hx[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[1,ip]
        hy[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[2,ip]
        hz[ip] = ( Vg_x[ip] + Vg_c[ip] ) * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip in 1:Npoints
        V_xc[ip] = V_x[ip] + V_c[ip] - 2.0*divh[ip]
        xc_calc.Vtau[ip,1] = Vtau_x[ip] + Vtau_c[ip]
    end

    return
end

