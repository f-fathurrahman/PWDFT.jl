"""
Initializes full ionic potential on a given `PWGrid`.
The structure factor `strf` must be calculated before.
Atomic coordinates information is thus contained in `strf`.
Charges of atomic nuclei `Znucls` for each species present in
atoms should also be provided.
"""
function init_V_coulomb_G( pw::PWGrid, strf::Array{ComplexF64,2}, Znucls::Array{Float64,1} )

    Nsp1 = size(strf)[2]
    Nsp2 = size(Znucls)[1]

    # Check for consistency between given no. pseudopots and no. of species
    # present in atoms.
    if Nsp1 != Nsp2
        error( @sprintf("Nsp1 /= Nsp2: %d %d\n", Nsp1, Nsp2) )
    end
    Nspecies = Nsp1

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    Vg = zeros(ComplexF64, Npoints)
    V  = zeros(Float64, Npoints)
    #
    for isp = 1:Nspecies
        prefactor = -4*pi*Znucls[isp]/CellVolume
        # Note that Vg[1] (for GVector zero) is always zero 0.0
        # To be sure let's put is here anyway.
        Vg[1] = 0.0 + im*0.0
        for ig = 2:Ng
            ip = idx_g2r[ig]
            Vg[ip] = prefactor/G2[ig]*strf[ig,isp]
        end
        V[:] = V[:] + real( G_to_R(pw, Vg) ) * Npoints
    end
    return V
end
