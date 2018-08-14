function calc_epsxc_VWN( Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    ptr = ccall( (:xc_func_alloc, LIBXC), Ptr{XCFuncType}, () )

    #
    # exchange part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Cint),
            ptr, 1, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, eps_x )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    # correlation part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
            ( Ptr{XCFuncType}, Cint, Cint),
            ptr, 7, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, eps_c )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    ccall( (:xc_func_free, LIBXC), Cvoid, (Ref{XCFuncType},), ptr )

    return eps_x + eps_c

end

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

    # Do transpose manually
    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)
    
    ptr = ccall( (:xc_func_alloc, LIBXC), Ptr{XCFuncType}, () )

    #
    # exchange part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Cint),
            ptr, 1, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe_tmp, eps_x )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    # correlation part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
            ( Ptr{XCFuncType}, Cint, Cint),
            ptr, 7, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe_tmp, eps_c )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    ccall( (:xc_func_free, LIBXC), Cvoid, (Ref{XCFuncType},), ptr )

    return eps_x + eps_c
end


function calc_Vxc_VWN( Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    v_x = zeros(Float64,Npoints)
    v_c = zeros(Float64,Npoints)

    ptr = ccall( (:xc_func_alloc, LIBXC), Ptr{XCFuncType}, () )

    #
    # exchange part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Cint),
            ptr, 1, Nspin)
    
    ccall( (:xc_lda_vxc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, v_x )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    # correlation part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
            ( Ptr{XCFuncType}, Cint, Cint),
            ptr, 7, Nspin)
    
    ccall( (:xc_lda_vxc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, v_c )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    ccall( (:xc_func_free, LIBXC), Cvoid, (Ref{XCFuncType},), ptr )

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
    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]
    
    ptr = ccall( (:xc_func_alloc, LIBXC), Ptr{XCFuncType}, () )

    #
    # exchange part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Cint),
            ptr, 1, Nspin)

    ccall( (:xc_lda_vxc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe_tmp, V_x )

    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    # correlation part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           ( Ptr{XCFuncType}, Cint, Cint),
           ptr, 7, Nspin)

    ccall( (:xc_lda_vxc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe_tmp, V_c )

    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    ccall( (:xc_func_free, LIBXC), Cvoid, (Ref{XCFuncType},), ptr )


    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = V_x[ip] + V_c[ip]
        Vxc[ipp,2] = V_x[ip+1] + V_c[ip+1]
    end
    return Vxc
end

