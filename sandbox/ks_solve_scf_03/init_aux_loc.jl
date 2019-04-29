include("Gshells.jl")

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

    Nelectrons = Ham.electrons.Nelectrons

    strf = calc_strfact( atoms, pw )

    G2_gshells, idx_gshells = init_Gshells( pw.gvec )
    ngl = length(G2_gshells)

    Vgl = zeros(Float64, ngl, Nspecies)

    for isp = 1:Nspecies
        psp = pspots[isp]
        for igl = 2:ngl
            Vgl[igl,isp] = eval_Vloc_G( psp, G2_gshells[igl] )/CellVolume
        end
    end

    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    for isp = 1:Nspecies
        for ig = 2:Ng
            ip = idx_g2r[ig]
            igl = idx_gshells[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl,isp]
        end
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    #

    E_alphat = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        #myfunc(r) = r^2 * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]/r )
        myfunc(r) = r^2 * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]*erf(r)/r )
        E_alphat = E_alphat + 4*pi*quadgk( myfunc, eps(), 10.0 )[1]*Zvals[isp]
    end

    E_alphat = E_alphat/CellVolume
    println("E_alphat = ", E_alphat)

    E_sisa = pi*Nelectrons^2/CellVolume
    println("E_sisa = ", E_sisa)


    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    Rhoe_gl = zeros(Float64,ngl)
    for igl = 1:ngl
        Rhoe_gl[igl] = exp( -0.25*G2_gshells[igl] )
    end
    for isp = 1:Nspecies
        for ig = 1:Ng
            ip = idx_g2r[ig]
            igl = idx_gshells[ig]
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] - Zvals[isp]*Rhoe_gl[igl]*strf[ig,isp]/CellVolume
        end
    end
    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )*Npoints


    println("integ Rhoe_aux = ", sum(Rhoe_aux)*CellVolume/Npoints)

    return V_Ps_loc, Rhoe_aux, E_alphat, E_sisa

end

