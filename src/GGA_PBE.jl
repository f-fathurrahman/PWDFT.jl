function calc_Vxc_PBE( pw::PWGrid, Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] +
                     gRhoe[3,ip]*gRhoe[3,ip]
    end

    Vxc = zeros( Float64, Npoints )
    Vgxc = zeros( Float64, Npoints )
    #
    ccall( (:calc_Vxc_PBE, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, gRhoe2, Vxc, Vgxc )

    h = zeros(Float64,3,Npoints)
    for ip = 1:Npoints
        h[1,ip] = Vgxc[ip] * gRhoe[1,ip]
        h[2,ip] = Vgxc[ip] * gRhoe[2,ip]
        h[3,ip] = Vgxc[ip] * gRhoe[3,ip]
    end

    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, h )

    for ip = 1:Npoints 
      Vxc[ip] = Vxc[ip] - 2.0*divh[ip]
    end

    return Vxc
end


function calc_epsxc_PBE( Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]
    epsxc = zeros( Float64, Npoints )

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] +
                     gRhoe[3,ip]*gRhoe[3,ip]
    end
    #
    ccall( (:calc_epsxc_PBE, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, gRhoe2, epsxc )
    #
    return epsxc
end



# Probably should be moved to PWGrid
function op_nabla( pw::PWGrid, Rhoe::Array{Float64,1} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    RhoeG = R_to_G(pw,Rhoe)[idx_g2r]

    ∇RhoeG_full = zeros(Complex128,3,Npoints)
    ∇Rhoe = zeros(Float64,3,Npoints)
    
    for ig = 1:Ng
        ip = idx_g2r[ig]
        ∇RhoeG_full[1,ip] = im*G[1,ig]*RhoeG[ig]
        ∇RhoeG_full[2,ip] = im*G[2,ig]*RhoeG[ig]
        ∇RhoeG_full[3,ip] = im*G[3,ig]*RhoeG[ig]
    end

    ∇Rhoe[1,:] = real(G_to_R(pw,∇RhoeG_full[1,:]))
    ∇Rhoe[2,:] = real(G_to_R(pw,∇RhoeG_full[2,:]))
    ∇Rhoe[3,:] = real(G_to_R(pw,∇RhoeG_full[3,:]))
    return ∇Rhoe

end


function op_nabla_dot( pw::PWGrid, h::Array{Float64,2} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    hG = zeros(Complex128,3,Ng)
    hG[1,:] = R_to_G( pw, h[1,:] )[idx_g2r]
    hG[2,:] = R_to_G( pw, h[2,:] )[idx_g2r]
    hG[3,:] = R_to_G( pw, h[3,:] )[idx_g2r]

    divhG_full = zeros(Complex128,Npoints)
    
    for ig = 1:Ng
        ip = idx_g2r[ig]
        divhG_full[ip] = im*( G[1,ig]*hG[ig] + G[2,ig]*hG[ig] + G[3,ig]*hG[ig] )
    end

    divh = real( G_to_R( pw, divhG_full ) )
    return divh

end
