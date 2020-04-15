function calc_forces_Ps_loc( Ham )
    F_Ps_loc = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_Ps_loc!( Ham, F_Ps_loc )
    return F_Ps_loc
end

function calc_forces_Ps_loc!( Ham::Hamiltonian, F_Ps_loc::Array{Float64,2} )
    calc_forces_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, F_Ps_loc )
    return
end

function calc_forces_Ps_loc!(
    atoms::Atoms, pw::PWGrid, pspots::Array{PsPot_GTH,1},
    Rhoe::Array{Float64,2},
    F_Ps_loc::Array{Float64,2}
)
    
    Ng = pw.gvec.Ng
    Natoms = atoms.Natoms

    G = pw.gvec.G
    G2 = pw.gvec.G2

    idx_g2r = pw.gvec.idx_g2r
    
    Npoints = prod(pw.Ns)
    
    Rhoe_tot = zeros(Float64,Npoints)
    Nspin = size(Rhoe)[2]
    for ispin = 1:Nspin
        for ip = 1:Npoints
            Rhoe_tot[ip] = Rhoe_tot[ip] + Rhoe[ip,ispin]
        end
    end

    RhoeG = R_to_G(pw, Rhoe_tot)/Npoints  # this normalization is required

    Ω = pw.CellVolume

    atpos = atoms.positions
    
    for ia = 1:Natoms
        
        isp = atoms.atm2species[ia]
        psp = pspots[isp]
        
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig = 2:Ng
            
            GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
            Sf = cos(GX) - im*sin(GX)

            Vg_ig = eval_Vloc_G( psp, G2[ig] )

            ip = idx_g2r[ig]
            F_Ps_loc[:,ia] = F_Ps_loc[:,ia] + real(im*G[:,ig]*Vg_ig*conj(RhoeG[ip])*Sf)
        end
    end

    return F_Ps_loc
end