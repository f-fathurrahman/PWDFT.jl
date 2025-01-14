# TODO: cumulate spin 1,2 as one function or dispatch by julia's type system

# FIXME: These functions calculate gradient of Rhoe on the fly
# Libxc.evaluate! has simpler interface.

function calc_epsxc_Vxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::AbstractVector{Float64},
    epsxc::AbstractVector{Float64},
    V_xc::AbstractVector{Float64}
)
    Npoints = size(Rhoe, 1)
    Nspin = 1

    # calculate gρ2
    gρ = op_nabla( pw, Rhoe )
    gρ2 = zeros( Float64, Npoints )
    for ip in 1:Npoints
        gρ2[ip] = gρ[1,ip]^2 + gρ[2,ip]^2 + gρ[3,ip]^2
    end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
 
    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
 
    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gρ2, eps_x, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gρ2, eps_c, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    @views epsxc[:] .= eps_x[:] .+ eps_c[:]

    @views V_xc[:] .= V_x[:] .+ V_c[:] # update V_xc (the output)

    # gradient correction
    Vg_xc = Vg_x + Vg_c
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip in 1:Npoints # using linear indexing
        hx[ip] = Vg_xc[ip] * gρ[1,ip]
        hy[ip] = Vg_xc[ip] * gρ[2,ip]
        hz[ip] = Vg_xc[ip] * gρ[3,ip]
    end
    # div ( vgrho * gρ )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip in 1:Npoints
        V_xc[ip] = V_xc[ip] - 2.0*divh[ip]
    end

    return
end


function calc_epsxc_Vxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::Array{Float64,2},
    epsxc::Vector{Float64},
    V_xc::Array{Float64,2}
)
    Nspin = size(Rhoe, 2)
    if Nspin == 1
        @views calc_epsxc_Vxc_PBE!( xc_calc, pw, Rhoe[:,1], epsxc[:], V_xc[:,1] )
        return
    end

    Npoints = size(Rhoe, 1)

    # calculate gρ2
    gρ_up = op_nabla( pw, Rhoe[:,1] )
    gρ_dn = op_nabla( pw, Rhoe[:,2] )
    gρ2 = zeros( Float64, 3*Npoints )
    ipp = 0
    for ip in 1:3:3*Npoints
        ipp = ipp + 1
        # up-up
        gρ2[ip]   = gρ_up[1,ipp]*gρ_up[1,ipp] + gρ_up[2,ipp]*gρ_up[2,ipp] + gρ_up[3,ipp]*gρ_up[3,ipp]
        # up-dn
        gρ2[ip+1] = gρ_up[1,ipp]*gρ_dn[1,ipp] + gρ_up[2,ipp]*gρ_dn[2,ipp] + gρ_up[3,ipp]*gρ_dn[3,ipp]
        # dn-dn
        gρ2[ip+2] = gρ_dn[1,ipp]*gρ_dn[1,ipp] + gρ_dn[2,ipp]*gρ_dn[2,ipp] + gρ_dn[3,ipp]*gρ_dn[3,ipp]
    end

    Rhoe_tmp = zeros(Float64, 2*Npoints)
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    eps_x = zeros(Float64, Npoints)
    eps_c = zeros(Float64, Npoints)

    V_x  = zeros(Float64, Npoints*2)
    V_c  = zeros(Float64, Npoints*2)

    Vg_x  = zeros(Float64, 3*Npoints)
    Vg_c  = zeros(Float64, 3*Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe_tmp, gρ2, eps_x, V_x, Vg_x)
    Libxc_xc_func_end(ptr)
    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe_tmp, gρ2, eps_c, V_c, Vg_c)
    Libxc_xc_func_end(ptr)
    #
    Libxc_xc_func_free(ptr)

    #
    # energy
    #
    epsxc[:] .= eps_x[:] + eps_c[:]

    #
    # The potential
    #
    # LDA contrib
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        V_xc[ipp,1] = V_x[ip] + V_c[ip]
        V_xc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end
    #
    # Now the gradient correction
    #
    Vg_xc = reshape(Vg_x + Vg_c, (3,Npoints))
    # Vg_xc[1,:] -> up-up
    # Vg_xc[2,:] -> up-dn
    # Vg_xc[3,:] -> dn-dn
    #
    divh = zeros(Float64, pw.Ns)
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    #
    # spin up
    for ip in 1:Npoints
        hx[ip] = 2*Vg_xc[1,ip]*gρ_up[1,ip] + Vg_xc[2,ip]*gρ_dn[1,ip]
        hy[ip] = 2*Vg_xc[1,ip]*gρ_up[2,ip] + Vg_xc[2,ip]*gρ_dn[2,ip]
        hz[ip] = 2*Vg_xc[1,ip]*gρ_up[3,ip] + Vg_xc[2,ip]*gρ_dn[3,ip]
    end
    #
    divh[:,:,:] .= op_nabla_dot( pw, hx, hy, hz )
    #
    # Add the correction (-divh) for spin up
    for ip in 1:Npoints
        V_xc[ip,1] = V_xc[ip,1] - divh[ip]
    end # using linear index for divh
    #
    # Spin down
    for ip in 1:Npoints
        hx[ip] = 2*Vg_xc[3,ip]*gρ_dn[1,ip] + Vg_xc[2,ip]*gρ_up[1,ip]
        hy[ip] = 2*Vg_xc[3,ip]*gρ_dn[2,ip] + Vg_xc[2,ip]*gρ_up[2,ip]
        hz[ip] = 2*Vg_xc[3,ip]*gρ_dn[3,ip] + Vg_xc[2,ip]*gρ_up[3,ip]
    end
    #
    divh[:,:,:] .= op_nabla_dot( pw, hx, hy, hz )
    # Add the correction (-divh) for spin down
    for ip in 1:Npoints
        V_xc[ip,2] = V_xc[ip,2] - divh[ip]
    end
    #
    return
end

function calc_epsxc_PBE(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe
)
    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints) # epsxc does not depend on spin
    calc_epsxc_PBE!(xc_calc, pw, Rhoe, epsxc)
    return epsxc
end

# In-place version
function calc_epsxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::AbstractVector{Float64},
    epsxc::Vector{Float64}
)

    Npoints = size(Rhoe, 1)
    Nspin = 1

    # calculate gρ2
    gρ = op_nabla( pw, Rhoe )
    gρ2 = zeros( Float64, Npoints )
    for ip in 1:Npoints
        gρ2[ip] = gρ[1,ip]^2 + gρ[2,ip]^2 + gρ[3,ip]^2
    end

    eps_x = zeros(Float64, Npoints)
    eps_c = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe, gρ2, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe, gρ2, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    epsxc[:] .= eps_x[:] .+ eps_c[:]

    return
end

# spin-polarized
function calc_epsxc_PBE(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::Array{Float64,2}
)
    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    calc_epsxc_PBE!(xc_calc, pw, Rhoe, epsxc)
    return epsxc
end

function calc_epsxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::Array{Float64,2},
    epsxc::Vector{Float64}
)
    Nspin = size(Rhoe, 2)
    if Nspin == 1
        @views calc_epsxc_PBE!( xc_calc, pw, Rhoe[:,1], epsxc )
        return
    end

    Npoints = size(Rhoe, 1)

    # calculate gρ2
    gρ_up = op_nabla( pw, Rhoe[:,1] )
    gρ_dn = op_nabla( pw, Rhoe[:,2] )
    gρ2 = zeros( Float64, 3*Npoints )
    ipp = 0
    for ip in 1:3:3*Npoints
        ipp = ipp + 1
        # up-up
        gρ2[ip]   = gρ_up[1,ipp]*gρ_up[1,ipp] + gρ_up[2,ipp]*gρ_up[2,ipp] + gρ_up[3,ipp]*gρ_up[3,ipp]
        # up-dn
        gρ2[ip+1] = gρ_up[1,ipp]*gρ_dn[1,ipp] + gρ_up[2,ipp]*gρ_dn[2,ipp] + gρ_up[3,ipp]*gρ_dn[3,ipp]
        # dn-dn
        gρ2[ip+2] = gρ_dn[1,ipp]*gρ_dn[1,ipp] + gρ_dn[2,ipp]*gρ_dn[2,ipp] + gρ_dn[3,ipp]*gρ_dn[3,ipp]
    end

    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe_tmp, gρ2, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc!(ptr, Npoints, Rhoe_tmp, gρ2, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    epsxc[:] .= eps_x[:] .+ eps_c[:]
    return
end

function calc_Vxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::AbstractVector{Float64},
    V_xc::AbstractVector{Float64}
)

    Npoints = size(Rhoe)[1]
    Nspin = 1

    # calculate gρ2
    gρ = op_nabla( pw, Rhoe )
    gρ2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gρ2[ip] = gρ[1,ip]*gρ[1,ip] + gρ[2,ip]*gρ[2,ip] + gρ[3,ip]*gρ[3,ip]
    end

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)

    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe, gρ2, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe, gρ2, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    # LDA part
    V_xc[:] .= V_x[:] .+ V_c[:]
    
    # gradient correction
    Vg_xc = Vg_x + Vg_c
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip in 1:Npoints
        hx[ip] = Vg_xc[ip] * gρ[1,ip]
        hy[ip] = Vg_xc[ip] * gρ[2,ip]
        hz[ip] = Vg_xc[ip] * gρ[3,ip]
    end
    # div ( vgrho * gρ )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip in 1:Npoints
        V_xc[ip] = V_xc[ip] - 2.0*divh[ip] # factor 2 of spin-degeneracy
    end

    return
end


# spinpol
function calc_Vxc_PBE(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::Array{Float64,2},
)
    Npoints = size(Rhoe, 1)
    Nspin = size(Rhoe, 2)
    @assert Nspin == 2
    V_xc = zeros(Float64, Npoints, Nspin)
    calc_Vxc_PBE!(xc_calc, pw, Rhoe, V_xc)
    return V_xc
end

function calc_Vxc_PBE!(
    xc_calc::LibxcXCCalculator,
    pw,
    Rhoe::Array{Float64,2},
    V_xc::Array{Float64,2}
)

    Nspin = size(Rhoe, 2)
    if Nspin == 1
        @views calc_Vxc_PBE!( xc_calc, pw, Rhoe[:,1], Vxc[:,1] )
        return
    end

    Npoints = size(Rhoe, 1)

    # calculate gρ2
    gρ_up = op_nabla( pw, Rhoe[:,1] ) # gρ = ∇⋅Rhoe
    gρ_dn = op_nabla( pw, Rhoe[:,2] )
    #
    gρ2 = zeros( Float64, 3*Npoints )
    ipp = 0
    for ip in 1:3:3*Npoints
        ipp = ipp + 1
        gρ2[ip]   = gρ_up[1,ipp]*gρ_up[1,ipp] + gρ_up[2,ipp]*gρ_up[2,ipp] + gρ_up[3,ipp]*gρ_up[3,ipp]
        gρ2[ip+1] = gρ_up[1,ipp]*gρ_dn[1,ipp] + gρ_up[2,ipp]*gρ_dn[2,ipp] + gρ_up[3,ipp]*gρ_dn[3,ipp]
        gρ2[ip+2] = gρ_dn[1,ipp]*gρ_dn[1,ipp] + gρ_dn[2,ipp]*gρ_dn[2,ipp] + gρ_dn[3,ipp]*gρ_dn[3,ipp]
    end

    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    V_x  = zeros(Float64, Npoints*2)
    V_c  = zeros(Float64, Npoints*2)
    Vg_x  = zeros(Float64, 3*Npoints)
    Vg_c  = zeros(Float64, 3*Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe_tmp, gρ2, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_vxc!(ptr, Npoints, Rhoe_tmp, gρ2, V_c, Vg_c)
    Libxc_xc_func_end(ptr)
    #
    Libxc_xc_func_free(ptr)
    #
    # The potential
    #
    # LDA contrib
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        V_xc[ipp,1] = V_x[ip] + V_c[ip]
        V_xc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end
    #
    # Now the gradient correction
    #
    Vg_xc = reshape(Vg_x + Vg_c, (3,Npoints))
    # Vg_xc[1,:] -> up-up
    # Vg_xc[2,:] -> up-dn
    # Vg_xc[3,:] -> dn-dn
    #
    divh = zeros(Float64, pw.Ns)
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    #
    # spin up
    for ip in 1:Npoints
        hx[ip] = 2*Vg_xc[1,ip]*gρ_up[1,ip] + Vg_xc[2,ip]*gρ_dn[1,ip]
        hy[ip] = 2*Vg_xc[1,ip]*gρ_up[2,ip] + Vg_xc[2,ip]*gρ_dn[2,ip]
        hz[ip] = 2*Vg_xc[1,ip]*gρ_up[3,ip] + Vg_xc[2,ip]*gρ_dn[3,ip]
    end
    #
    divh[:,:,:] .= op_nabla_dot( pw, hx, hy, hz )
    #
    # Add the correction (-divh) for spin up
    for ip in 1:Npoints
        V_xc[ip,1] = V_xc[ip,1] - divh[ip]
    end
    #
    # Spin down
    for ip in 1:Npoints
        hx[ip] = 2*Vg_xc[3,ip]*gρ_dn[1,ip] + Vg_xc[2,ip]*gρ_up[1,ip]
        hy[ip] = 2*Vg_xc[3,ip]*gρ_dn[2,ip] + Vg_xc[2,ip]*gρ_up[2,ip]
        hz[ip] = 2*Vg_xc[3,ip]*gρ_dn[3,ip] + Vg_xc[2,ip]*gρ_up[3,ip]
    end
    #
    divh[:,:,:] .= op_nabla_dot( pw, hx, hy, hz )
    # Add the correction (-divh) for spin down
    for ip in 1:Npoints
        V_xc[ip,2] = V_xc[ip,2] - divh[ip]
    end
    #
    return
end
