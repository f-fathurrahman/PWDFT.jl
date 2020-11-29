function calc_forces_Ps_loc!( Ham::HamiltonianGamma, F_Ps_loc::Array{Float64,2} )
    calc_forces_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, F_Ps_loc )
    return
end

function calc_forces_Ps_loc!(
    atoms::Atoms, pw::PWGridGamma, pspots::Array{PsPot_GTH,1},
    Rhoe::Array{Float64,2},
    F_Ps_loc::Array{Float64,2}
)
    
    Ng = pw.gvec.Ng
    Natoms = atoms.Natoms

    G = pw.gvec.G
    G2 = pw.gvec.G2

    idx_g2r = pw.gvec.idx_g2r
    
    Npoints = prod(pw.Ns)
    # Prepare to transform Rhoe to RhoeG
    RhoeG = zeros(ComplexF64,pw.Ns)
    Nspin = size(Rhoe)[2]
    for ispin = 1:Nspin
        for ip = 1:Npoints
            RhoeG[ip] = RhoeG[ip] + Rhoe[ip,ispin]
        end
    end
    R_to_G!(pw, RhoeG)
    lmul!(1/Npoints, RhoeG) # this normalization is required

    atpos = atoms.positions
    for ia = 1:Natoms
        isp = atoms.atm2species[ia]
        psp = pspots[isp]
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig = 2:Ng
            # Structure factor, recalculate here
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf = cos(GX) - im*sin(GX)
            # Ps_loc, recalculate here
            Vg_ig = eval_Vloc_G(psp, G2[ig])
            ip = idx_g2r[ig]
            for i in 1:3
                F_Ps_loc[i,ia] = F_Ps_loc[i,ia] + 2*real(im*G[i,ig]*Vg_ig*conj(RhoeG[ip])*Sf)
            end
        end
    end
    return
end