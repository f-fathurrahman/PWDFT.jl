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
    # XXX Use sum or reduce
    for ispin in 1:Nspin, ip in 1:Npoints
        Rhoe_tot[ip] = Rhoe_tot[ip] + Rhoe[ip,ispin]
    end
    #
    RhoeG = R_to_G(pw, Rhoe_tot)/Npoints  # this normalization is required
    X = atoms.positions
    #
    for ia in 1:Natoms
        #   
        isp = atoms.atm2species[ia]
        psp = pspots[isp]
        #
        # G=0 contribution is zero, G[:,ig=0] = [0, 0, 0]
        for ig in 2:Ng
            #
            GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
            Sf = cos(GX) - im*sin(GX)
            #
            Vg_ig = eval_Vloc_G( psp, G2[ig] )
            #
            ip = idx_g2r[ig]
            @views F_Ps_loc[:,ia] += real(im*G[:,ig]*Vg_ig*conj(RhoeG[ip])*Sf)
        end
    end
    #
    return
end