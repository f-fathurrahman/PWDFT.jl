function init_V_coulomb_G( pw::PWGrid, strf::Array{Complex128,2}, Znucls::Array{Float64} )

    Nsp1 = size(strf)[2]
    Nsp2 = size(Znucls)[1]
    #
    # This is the restriction for this function
    #
    if Nsp1 != Nsp2
        @printf("ERROR: Nsp1 /= Nsp2: %d %d\n", Nsp1, Nsp2)
        exit()
    end
    Nspecies = Nsp1

    Npoints = prod(pw.Ns)
    Ω = pw.Ω

    Vg = zeros(Complex128, Npoints)
    G2 = pw.gvec.G2
    V  = zeros(Float64, Npoints)
    #
    for isp = 1:Nspecies
        prefactor = -4*pi*Znucls[isp]/Ω
        #
        # V[ig=0] = 0.0
        #
        for ig = 2:Npoints
            Vg[ig] = prefactor/G2[ig]*strf[ig,isp]
        end
        V[:] = V[:] + real( G_to_R(pw.Ns, Vg) ) * Npoints
    end
    return V
end
