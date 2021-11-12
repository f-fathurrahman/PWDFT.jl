using Libxc

# TODO cumulate spin 1,2 as one function or dispatch by julia's type system

function calc_epsxc_PBE( xc_calc::LibxcXCCalculator, pw, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe, gRhoe2, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_epsxc_PBE( xc_calc::LibxcXCCalculator, pw, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_epsxc_PBE( xc_calc, pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe2 = zeros( Float64, 3*Npoints )
    ipp = 0
    for ip = 1:3:3*Npoints
        ipp = ipp + 1
        gRhoe2[ip]   = gRhoe_up[1,ipp]*gRhoe_up[1,ipp] + gRhoe_up[2,ipp]*gRhoe_up[2,ipp] + gRhoe_up[3,ipp]*gRhoe_up[3,ipp]
        gRhoe2[ip+1] = gRhoe_up[1,ipp]*gRhoe_dn[1,ipp] + gRhoe_up[2,ipp]*gRhoe_dn[2,ipp] + gRhoe_up[3,ipp]*gRhoe_dn[3,ipp]
        gRhoe2[ip+2] = gRhoe_dn[1,ipp]*gRhoe_dn[1,ipp] + gRhoe_dn[2,ipp]*gRhoe_dn[2,ipp] + gRhoe_dn[3,ipp]*gRhoe_dn[3,ipp]
    end

    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe_tmp, gRhoe2, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe_tmp, gRhoe2, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    return eps_x + eps_c

end

function calc_Vxc_PBE( xc_calc::LibxcXCCalculator, pw, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
    V_xc = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
    Vg_xc = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe, gRhoe2, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    V_xc = V_x + V_c
    Vg_xc = Vg_x + Vg_c

    # gradient correction
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip = 1:Npoints
        hx[ip] = Vg_xc[ip] * gRhoe[1,ip]
        hy[ip] = Vg_xc[ip] * gRhoe[2,ip]
        hz[ip] = Vg_xc[ip] * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip = 1:Npoints
        V_xc[ip] = V_xc[ip] - 2.0*divh[ip]
    end


    return V_xc

end

function calc_Vxc_PBE( xc_calc::LibxcXCCalculator, pw, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_Vxc_PBE( xc_calc, pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] ) # gRhoe = ∇⋅Rhoe
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    #
    gRhoe2 = zeros( Float64, 3*Npoints )
    ipp = 0
    for ip = 1:3:3*Npoints
        ipp = ipp + 1
        gRhoe2[ip]   = gRhoe_up[1,ipp]*gRhoe_up[1,ipp] + gRhoe_up[2,ipp]*gRhoe_up[2,ipp] + gRhoe_up[3,ipp]*gRhoe_up[3,ipp]
        gRhoe2[ip+1] = gRhoe_up[1,ipp]*gRhoe_dn[1,ipp] + gRhoe_up[2,ipp]*gRhoe_dn[2,ipp] + gRhoe_up[3,ipp]*gRhoe_dn[3,ipp]
        gRhoe2[ip+2] = gRhoe_dn[1,ipp]*gRhoe_dn[1,ipp] + gRhoe_dn[2,ipp]*gRhoe_dn[2,ipp] + gRhoe_dn[3,ipp]*gRhoe_dn[3,ipp]
    end

    V_xc = zeros(Float64, Npoints, 2)
    V_x  = zeros(Float64, Npoints*2)
    V_c  = zeros(Float64, Npoints*2)

    Vg_xc = zeros(Float64, 3, Npoints)
    Vg_x  = zeros(Float64, 3*Npoints)
    Vg_c  = zeros(Float64, 3*Npoints)

    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end


    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe_tmp, gRhoe2, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe_tmp, gRhoe2, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        V_xc[ipp,1] = V_x[ip] + V_c[ip]
        V_xc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end

    Vg_xc = reshape(Vg_x + Vg_c, (3,Npoints))

    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    divh = zeros(Float64,Npoints)

    #
    # spin up
    #
    for ip = 1:Npoints
        hx[ip] = 2*Vg_xc[1,ip]*gRhoe_up[1,ip] + Vg_xc[2,ip]*gRhoe_dn[1,ip]
        hy[ip] = 2*Vg_xc[1,ip]*gRhoe_up[2,ip] + Vg_xc[2,ip]*gRhoe_dn[2,ip]
        hz[ip] = 2*Vg_xc[1,ip]*gRhoe_up[3,ip] + Vg_xc[2,ip]*gRhoe_dn[3,ip]
    end

    divh = op_nabla_dot( pw, hx, hy, hz )

    # spin up
    for ip = 1:Npoints
        V_xc[ip,1] = V_xc[ip,1] - divh[ip]
    end

    #
    # Spin down
    #
    for ip = 1:Npoints
        hx[ip] = 2*Vg_xc[3,ip]*gRhoe_dn[1,ip] + Vg_xc[2,ip]*gRhoe_up[1,ip]
        hy[ip] = 2*Vg_xc[3,ip]*gRhoe_dn[2,ip] + Vg_xc[2,ip]*gRhoe_up[2,ip]
        hz[ip] = 2*Vg_xc[3,ip]*gRhoe_dn[3,ip] + Vg_xc[2,ip]*gRhoe_up[3,ip]
    end

    #
    divh = op_nabla_dot( pw, hx, hy, hz )
    # spin down
    for ip = 1:Npoints
        V_xc[ip,2] = V_xc[ip,2] - divh[ip]
    end

    return V_xc

end
