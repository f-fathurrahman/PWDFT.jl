function calc_epsxc_LDA( Rhoe::Array{Float64,1}, xc_id )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    epsxc = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()

    Libxc.xc_func_init(ptr, xc_id, Nspin)
    Libxc.xc_lda_exc!(ptr, Npoints, Rhoe, epsxc)
    Libxc.xc_func_end(ptr)

    Libxc.xc_func_free(ptr)

    return epsxc
end


#=
"""
Calculate XC energy per particle using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_epsxc_VWN( Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    Npoints = size(Rhoe)[1]

    if Nspin == 1
        return calc_epsxc_VWN( Rhoe[:,1] )
    end

    # Do the transpose manually
    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end


    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.LDA_X, Nspin)
    Libxc.xc_lda_exc!(ptr, Npoints, Rhoe_tmp, eps_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.LDA_C_VWN, Nspin)
    Libxc.xc_lda_exc!(ptr, Npoints, Rhoe_tmp, eps_c)
    Libxc.xc_func_end(ptr)

    #
    Libxc.xc_func_free(ptr)

    return eps_x + eps_c
end


function calc_Vxc_VWN( Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    v_x = zeros(Float64,Npoints)
    v_c = zeros(Float64,Npoints)

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.LDA_X, Nspin)
    Libxc.xc_lda_vxc!(ptr, Npoints, Rhoe, v_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.LDA_C_VWN, Nspin)
    Libxc.xc_lda_vxc!(ptr, Npoints, Rhoe, v_c)
    Libxc.xc_func_end(ptr)

    #
    Libxc.xc_func_free(ptr)

    return v_x + v_c

end

"""
Calculate XC potential using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_Vxc_VWN( Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_Vxc_VWN( Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    Vxc = zeros( Float64, Npoints, 2 )
    V_x = zeros( Float64, 2*Npoints )
    V_c = zeros( Float64, 2*Npoints )

    # This is the transposed version of Rhoe, use copy
    Rhoe_tmp = zeros(2*Npoints)
    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Rhoe_tmp[ip] = Rhoe[ipp,1]
        Rhoe_tmp[ip+1] = Rhoe[ipp,2]
    end

    ptr = Libxc.xc_func_alloc()
    # exchange part
    Libxc.xc_func_init(ptr, Libxc.LDA_X, Nspin)
    Libxc.xc_lda_vxc!(ptr, Npoints, Rhoe_tmp, V_x)
    Libxc.xc_func_end(ptr)

    #
    # correlation part
    Libxc.xc_func_init(ptr, Libxc.LDA_C_VWN, Nspin)
    Libxc.xc_lda_vxc!(ptr, Npoints, Rhoe_tmp, V_c)
    Libxc.xc_func_end(ptr)

    #
    Libxc.xc_func_free(ptr)

    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = V_x[ip] + V_c[ip]
        Vxc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end
    return Vxc
end
=#
