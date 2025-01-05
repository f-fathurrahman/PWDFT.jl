# XXX where is this used?
function calc_epsxc_Vxc_VWN(
    xc_calc::LibxcXCCalculator,
    Rhoe::Array{Float64,2}
)
    Nspin = size(Rhoe, 2)
    @assert Nspin == 1
    if Nspin == 1
        return calc_epsxc_Vxc_VWN( xc_calc, Rhoe[:,1] )
    end
end

# In-place version, both epsxc and Vxc, spinpol
function calc_epsxc_Vxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::Array{Float64,2},
    epsxc::Vector{Float64},
    Vxc::Array{Float64,2}
)

    Nspin = size(Rhoe, 2)
    Npoints = size(Rhoe, 1)

    if Nspin == 1
        # Don't forget @views because Vxc will be modified
        @views calc_epsxc_Vxc_VWN!(xc_calc, Rhoe[:,1], epsxc, Vxc[:,1])
        return
    end

    # Do the transpose manually
    Rhoe_tmp = zeros(Float64, 2*Npoints)
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    eps_x = zeros(Float64, Npoints)
    eps_c = zeros(Float64, Npoints)

    V_x = zeros(Float64, 2*Npoints)
    V_c = zeros(Float64, 2*Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc_vxc!(ptr, Npoints, Rhoe_tmp, eps_x, V_x)
    Libxc_xc_func_end(ptr)
    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc_vxc!(ptr, Npoints, Rhoe_tmp, eps_c, V_c)
    Libxc_xc_func_end(ptr)
    #
    Libxc_xc_func_free(ptr)

    # Set outputs
    epsxc[:] .= eps_x[:] .+ eps_c[:]
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = V_x[ip] + V_c[ip]
        Vxc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end

    return
end


function calc_epsxc_Vxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::AbstractVector{Float64},
    epsxc::AbstractVector{Float64},
    Vxc::AbstractVector{Float64}
)

    Npoints = size(Rhoe, 1)
    Nspin = 1
    eps_x = zeros(Float64, Npoints)
    eps_c = zeros(Float64, Npoints)
    v_x = zeros(Float64, Npoints)
    v_c = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)  # LDA_X
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc_vxc!(ptr, Npoints, Rhoe, eps_x, v_x)
    Libxc_xc_func_end(ptr)
    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin) # LDA_C_VWN
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc_vxc!(ptr, Npoints, Rhoe, eps_c, v_c)
    Libxc_xc_func_end(ptr)
    #
    Libxc_xc_func_free(ptr)

    epsxc[:] .= eps_x[:] .+ eps_c[:]
    Vxc[:] .= v_x[:] .+ v_c[:]

    return
end


function calc_epsxc_Vxc_VWN(
    xc_calc::LibxcXCCalculator,
    Rhoe::AbstractVector{Float64}
)
    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    Vxc = zeros(Float64, Npoints)
    calc_epsxc_Vxc_VWN!( xc_calc, Rhoe, epsxc, Vxc )
    return epsxc, Vxc
end

#
# epsxc only
#
function calc_epsxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::AbstractVector{Float64},
    epsxc::Vector{Float64}
)

    Npoints = size(Rhoe)[1]
    Nspin = 1
    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)  # LDA_X
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe, eps_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin) # LDA_C_VWN
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    @views epsxc[:] = eps_x[:] + eps_c[:]
    return
end

function calc_epsxc_VWN(
    xc_calc::LibxcXCCalculator,
    Rhoe
)
    epsxc = zeros(Float64, size(Rhoe)) # always Vector{Float64}
    calc_epsxc_VWN!(xc_calc, Rhoe, epsxc)
    return epsxc
end



"""
Calculate XC energy per particle using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_epsxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::Array{Float64,2},
    epsxc::Vector{Float64}
)

    Nspin = size(Rhoe, 2)
    if Nspin == 1
        calc_epsxc_VWN!( xc_calc, Rhoe[:,1], epsxc )
        return
    end

    Npoints = size(Rhoe, 1)

    # Do the transpose manually
    Rhoe_tmp = zeros(Float64, 2*Npoints)
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    eps_x = zeros(Float64, Npoints)
    eps_c = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe_tmp, eps_x)
    Libxc_xc_func_end(ptr)
    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe_tmp, eps_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    epsxc[:] .= eps_x[:] .+ eps_c[:]
end


function calc_Vxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::AbstractVector{Float64},
    Vxc::AbstractVector{Float64}
)

    Npoints = size(Rhoe)[1]
    Nspin = 1
    v_x = zeros(Float64,Npoints)
    v_c = zeros(Float64,Npoints)

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10) # set in QE
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe, v_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe, v_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    @views Vxc[:] = v_x[:] + v_c[:]
    return
end


function calc_Vxc_VWN(
    xc_calc::LibxcXCCalculator,
    Rhoe::AbstractVector{Float64}
)
    Vxc = zeros(Float64, size(Rhoe))
    calc_Vxc_VWN!(xc_calc, Rhoe, Vxc)
    return Vxc
end




"""
Calculate XC potential using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_Vxc_VWN!(
    xc_calc::LibxcXCCalculator,
    Rhoe::Array{Float64,2},
    Vxc::Array{Float64,2}
)
    Nspin = size(Rhoe, 2)
    Npoints = size(Rhoe, 1)

    if Nspin == 1
        @views calc_Vxc_VWN!( xc_calc, Rhoe[:,1], Vxc[:,1] )
        return
    end

    V_x = zeros(Float64, 2*Npoints)
    V_c = zeros(Float64, 2*Npoints)

    # This is the transposed version of Rhoe, use copy
    Rhoe_tmp = zeros(Float64, 2*Npoints)
    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, xc_calc.x_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe_tmp, V_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, xc_calc.c_id, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_lda_vxc!(ptr, Npoints, Rhoe_tmp, V_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    ipp = 0
    for ip in 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = V_x[ip] + V_c[ip]
        Vxc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end
    return
end
