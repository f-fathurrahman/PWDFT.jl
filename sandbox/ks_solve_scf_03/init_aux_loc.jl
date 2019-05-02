function init_aux_loc( Ham )

    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    Zvals = atoms.Zvals

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2

    Nelectrons = Ham.electrons.Nelectrons

    gcut = 2.0
    ebsl = 1e-8
    glast2 = gcut*gcut
    gexp = -log(ebsl)    
    η = sqrt(glast2/gexp)/2

    strf = calc_strfact( atoms, pw )

    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    for isp = 1:Nspecies
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] - Zvals[isp]*exp(-0.125*G2[ig]/η^2) * strf[ig,isp] / CellVolume
        end
    end
    Rhoe_aux = real( G_to_R(pw, Rhoe_aux_G) )*Npoints
    println("integ Rhoe_aux = ", sum(Rhoe_aux)*CellVolume/Npoints)

    u_Ps_loc = zeros(Float64,Npoints)
    Vg = zeros(ComplexF64,Npoints)
    for isp = 1:Nspecies
        psp = pspots[isp]
        for ig = 2:Ng
            ip = idx_g2r[ig]
            Vg[ip] = ( eval_Vloc_G(psp, G2[ig]) + 4*pi*Zvals[isp]*exp(-0.125*G2[ig]/η^2)/G2[ig] ) * strf[ig,isp] / CellVolume
        end
        u_Ps_loc[:] = u_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end
    println("integ u_Ps_loc = ", sum(u_Ps_loc)*CellVolume/Npoints)

    return Rhoe_aux, u_Ps_loc

end

